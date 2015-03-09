"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    © 2014-2015. All rights reserved.
"""
"""
This Blender script that generates cross-sectional data from a skeletonization
representation in an Amiramesh text file, and a Blender project containing the
mesh for the corresponding skeletonized object.

NOTES:  It is important that the Blender project file *.blend contains the named object; also,
            The object must be in the correct position (typically this is the object whose
            mesh was used to create the skeletonization).
        The coordinate systems differ between Blender and Avizo (the script accounts for this).
        Blender doesn't release deleted meshes, so it is better not to chunk multiple nodes (memory usage drastically grows)
        Blender will fail if it doesn't get all the cores it expects, use -t 1, and don't run more copies of Blender than real cores.
        For now, run a one node test run on the cell to get the total number of nodes (part of the output file name).
            If there are N nodes, then use `echo $(seq 0 1 $(( N-1 )))` to iterate over the node chunks.
        Each invocation of the script creates a tab-delimited CSV file in the same directory as the skeleton *.am file
            The file is name after the chunk of segment cross-sectional data it contains
            Combine these files together after running all the scripts.
        Below are recipes for running this script as a single invocation, and in parallel.
        There are also some testing functions that can be run interactively (below), to verify that the script is working
            I.e., The cross sectional cuts have the correct cutting planes and cross-sectional faces.
            See these functions for more documentation.
        For interactive development from the Blender Python Console, enable reloading of the skeleton_annotate module as follows:
            import imp
            imp.reload(skeleton_annotate)

# single invocation of script
blender -t 1 -b <path>/<filename>.blend -P <path>/skeleton_annotate.py -- "<object_name>" /<path>/<skeleton_filename>.am <start_node_idx> <number_of_nodes>

# parallel invocation of script
echo $(seq 0 1 $(( <total_node_count> - 1)) ) | xargs -d " " -n 1 -P <num_cores> -I{} sh -c 'blender -t 1 -b <path>/<filename>.blend -P <path>/skeleton_annotate.py -- "<object_name>" /<path>/<skeleton_filename>.am $(( {} )) 1'

# combine *.csv data files into a single file (file is named *.dat so it isn't re-read recursively; rename to *.csv file afterwards)
head -n +1 "`ls *.csv | head -1`" > cross_section.csv.dat
for i in *.csv ; do tail -n +2 "$i" >> cross_section.csv.dat ; done
"""

import os
import sys
import csv
import json
import re
import logging
import operator
import math
import functools
from collections import defaultdict

import bmesh
import bpy
from bpy.props import *
import mathutils
import addon_utils

try:
    import skeletonizer
except ImportError:
    sys.path.append(os.path.abspath(os.path.split(__file__)[0]))

from skeletonizer.amiramesh import *

# TODO: Check for available object_cross_section addon
'''
print(bpy.context.user_preferences.addons.keys())
for mod in addon_utils.modules():
    print(mod)
    #print(mod.bl_info.get('version', (-1, -1, -1)))
if 'object_cross_section' not in bpy.context.user_preferences.addons.keys():
    addon_utils.enable('object_cross_section', default_set=True, persistent=False, handle_error=None)

#from object_cross_section import *
'''


def get_paths(arg):
    skel_path = os.path.abspath(os.path.dirname(arg))
    if arg[-3:] == '.am':
        skel_name = os.path.basename(arg[:-3])
    else:
        skel_name = os.path.basename(arg)

    skel_am_file = os.path.join(skel_path, skel_name + '.am')
    skel_json_file = os.path.join(skel_path, skel_name + '.annotations.json')

    return (skel_path, skel_name, skel_am_file, skel_json_file)


def calc_rotation(v):
    '''
    Calculates the XYZ rotation from a normal vector (tuple)
    :param v: The normal vector
    :return: The Euler rotation for an object for face normal vector v
    '''
    return mathutils.Vector((0,0,1)).rotation_difference(mathutils.Vector(v)).to_euler()


def swizzle_coordinates(p):
    '''
    XYZ coordinates in Avizo become -XZY in Blender, and visa-versa.
    Converts point in one coordinate system to the other.
    :param p: Input position 3-tuple.
    :return: Position converted into new coordinate system
    '''
    return (-p[0], p[2], p[1])



def generate_node_cross_section_data(obj_name, npos, nnorm):
    '''
    Generate cross-sectional data for object at a node position
    :param obj_name: Object name (string)
    :param npos: Node position vector (Vector), Blender coordinates
    :param nnorm: Node normal vector (Vector), Blender coordinates
    :return: Cross-sectional data tuple
    '''
    plane_ob_name = 'Plane'
    pos_ob_name = 'Partofsection'
    assert (plane_ob_name not in bpy.data.objects), "Expected special name '%s' to be unused" % (plane_ob_name)
    assert (pos_ob_name not in bpy.data.objects), "Expected special name '%s' to be unused" % (pos_ob_name)
    assert (obj_name in bpy.data.objects), "Expected object with name: '%s'" % (obj_name)

    bpy.ops.mesh.primitive_plane_add(location=npos, rotation=calc_rotation(nnorm))

    cell_ob = bpy.data.objects[obj_name]
    plane_ob = bpy.data.objects[plane_ob_name]

    cell_ob.hide = False
    bpy.ops.object.select_all(False)
    cell_ob.select = True
    plane_ob.select = True

    bpy.ops.object.cross_section()

    cell_ob.hide = True

    ob = bpy.data.objects[pos_ob_name]
    bpy.ops.object.select_all(False)

    plane_ob.select = True
    bpy.ops.object.delete()

    bm = bmesh.new()
    me = bm.from_mesh(ob.data)
    cx_faces = [(f, (npos - (ob.location + f.calc_center_median_weighted())).length_squared) for f in bm.faces]
    cx = functools.reduce(lambda a, b: b if not a else b if b[1] < a[1] else a, cx_faces, None)

    cx_area = cx[0].calc_area()
    cx_perim = cx[0].calc_perimeter()

    bm.free()

    ob.select = True
    bpy.ops.object.delete()

    return {'area':cx_area, 'perimeter':cx_perim}


def generate_cross_sections(obj_name, skel_am_file, skel_json_file, segment_range,
                            out_path):
    '''
    Generate cross section annotation data for object specified by skeleton in segment range.
    Creates a *.csv file containing per-segment point cross-sectional data
    :param obj_name: The Blender object name containing the corresponding mesh.
    :param skel_am_file: The Amiramesh file containing the skeleton
    :param skel_json_file: The annotation file for the corresponding skeleton
    :param segment_range: A tuple with the [start, end) range.
    :param out_path: Directory path for output file
    :return: A dictionary containing the cross-sectional data mapped to the
             original 3D segment point position in the *.am file.
    '''
    cxs = {}
    reader = AmirameshReader()
    with open(skel_am_file, 'r') as f:
        skel = reader.parse(f)
    with open(skel_json_file, 'r') as f:
        data = json.load(f)

    r = (max(0, segment_range[0]), min(len(skel.segments), segment_range[1]))

    csvfilename = obj_name+'-cross_section_data-range-'+str(r[0])+'-'+str(r[1]-1)+'-of-'+str(len(skel.segments))
    csv_file = os.path.join(out_path, csvfilename + '.csv')

    with open(csv_file, 'w', newline='') as f:
        csvheader = ['am_position','segment_idx','pnt_idx',
                     'area', 'perimeter',
                     'estimated_diameter', 'estimated_area', 'estimated_perimeter',
                     'blender_position', 'blender_normal']
        writer = csv.DictWriter(f, fieldnames=csvheader, delimiter='\t', quotechar='|')
        writer.writeheader()

        for idx in range(r[0], r[1]):
            s = skel.segments[idx]
            am_pts = [i for i in s.points]
            if len(am_pts) >= 2:
                prev_pnt = mathutils.Vector(swizzle_coordinates(am_pts[1].position()))
                p_idx = 0
                for p in am_pts:
                    am_ppos = p.position()
                    ppos = swizzle_coordinates(am_ppos)
                    pnt = mathutils.Vector(ppos)
                    pnorm = pnt - prev_pnt

                    cx_data = generate_node_cross_section_data(obj_name, pnt, pnorm)
                    n_data = {'segment_idx':idx, "pnt_idx":p_idx,
                              'am_position':am_ppos, 'blender_position':ppos, 'blender_normal':pnorm,
                              'estimated_diameter':p.diameter,
                              'estimated_area':math.pi * ((p.diameter / 2.0)**2),
                              'estimated_perimeter':math.pi * p.diameter}
                    cx_data.update(n_data)
                    cxs[am_ppos] = cx_data

                    writer.writerow(cx_data)
                    f.flush()
                    print("processed data:%s" % (cx_data))

                    p_idx += 1

    return cxs


def debug_cut_planes(obj_name, skeleton, segment_range):
    '''
    Creates cut planes for major nodes in the given skeletonization.
    '''

    r = (max(0, segment_range[0]), min(len(skeleton.segments), segment_range[1]))

    for idx in range(r[0], r[1]):
        s = skeleton.segments[idx]
        am_pts = [i for i in s.points]
        if len(am_pts) >= 2:
            prev_pnt = mathutils.Vector(swizzle_coordinates(am_pts[1].position()))
            p_idx = 0
            for p in am_pts:
                am_ppos = p.position()
                ppos = swizzle_coordinates(am_ppos)
                pnt = mathutils.Vector(ppos)
                pnorm = pnt - prev_pnt

                bpy.ops.mesh.primitive_plane_add(location=pnt, rotation=calc_rotation(pnorm))

                p_idx += 1
                if r[1] - r[0] > 1:
                    break   # for a single segment, create all cutting planes; otherwise, only start pnt

def debug_cut_faces(obj_name, skeleton, segment_range):
    '''
    Creates cut faces for major nodes in the given skeletonization.
    '''
    plane_ob_name = 'Plane'
    pos_ob_name = 'Partofsection'

    r = (max(0, segment_range[0]), min(len(skeleton.segments), segment_range[1]))

    for idx in range(r[0], r[1]):
        s = skeleton.segments[idx]
        am_pts = [i for i in s.points]
        if len(am_pts) >= 2:
            prev_pnt = mathutils.Vector(swizzle_coordinates(am_pts[1].position()))
            p_idx = 0
            for p in am_pts:
                am_ppos = p.position()
                ppos = swizzle_coordinates(am_ppos)
                pnt = mathutils.Vector(ppos)
                pnorm = pnt - prev_pnt

                bpy.ops.mesh.primitive_plane_add(location=pnt, rotation=calc_rotation(pnorm))

                cell_ob = bpy.data.objects[obj_name]
                plane_ob = bpy.data.objects[plane_ob_name]


                cell_ob.hide = False
                bpy.ops.object.select_all(False)
                cell_ob.select = True
                plane_ob.select = True

                bpy.ops.object.cross_section()

                cell_ob.hide = True

                ob = bpy.data.objects[pos_ob_name]
                bpy.ops.object.select_all(False)

                bm = bmesh.new()
                bm.from_mesh(ob.data)
                cx_faces = [(f, (pnt - (ob.location + f.calc_center_median_weighted())).length_squared) for f in bm.faces]
                cx = functools.reduce(lambda a, b: b if not a else b if b[1] < a[1] else a, cx_faces, None)
                cx_others = [f for f in bm.faces if not f is cx[0]]
                bmesh.ops.delete(bm, geom=cx_others, context=5) # context 5 is DEL_FACES
                bm.to_mesh(ob.data)
                bm.free()

                ob.name = 'Testslice'
                ob.select = False
                plane_ob.name = 'Testplane'
                plane_ob.select = False

                p_idx += 1
                if r[1] - r[0] > 1:
                    break   # for a single segment, create all cutting planes; otherwise, only start pnt



def test_cut_planes(segment_range = (0,1000)):
    '''
    A test function to show all the cutting planes used to create the cross-sections.
    After running this function the scene will contain 'Plane*' objects centered on
    skeleton nodes and oriented as normal cuts.  Visually verify the location and orientation.

    To run interactively, start blender with the project:
        blender <project_file>.blend
    Open the Python console and type:
        import sys, os
        sys.path.append(os.path.abspath('/var/remote/projects/epfl/development/staging/skeletonizer/Source'))
        import skeleton_annotate
    Run with:
        skeleton_annotate.test_cut_planes()

    NOTE: Currently hardcoded to 'Astrocyte 2' in KB-E0010.blend stack
    :param segment_range: A (start, end) tuple specifying which segment nodes to process.
                          A range of (n,n+1) will create cut planes for all the points in the segment,
                          otherwise the cuts are created for each start segment node.
    '''

    dpath = '/var/remote/projects/epfl/data/KB-E0010/astrocyte2'
    skel_am_file = dpath + '/GeometrySurface.Smt.SptGraph_OK.am'
    skel_json_file = dpath+'/GeometrySurface.Smt.SptGraph_OK.annotations.json'
    obj_name = 'Astrocyte 2'

    reader = AmirameshReader()
    with open(skel_am_file, 'r') as f:
        skeleton = reader.parse(f)
    with open(skel_json_file, 'r') as f:
        data = json.load(f)

    debug_cut_planes(obj_name, skeleton, segment_range)

def test_cut_faces(segment_range = (0,10)):
    '''
    A test function to show all the cutting planes and faces used to create the cross-sections.
    After running this function the scene will contain 'Testplane*' and 'Testslice*' objects centered on
    skeleton nodes and oriented as normal cuts and faces.
    Visually verify the location, orientation, and correctness of the selected face.
    'Testplane###' and 'Testslice###' should pair for given '###'.

    To run interactively, start blender with the project:
        blender <project_file>.blend
    Open the Python console and type:
        import sys, os
        sys.path.append(os.path.abspath('/var/remote/projects/epfl/development/staging/skeletonizer/Source'))
        import skeleton_annotate
    Run with:
        skeleton_annotate.test_cut_faces()

    NOTE: Currently hardcoded to 'Astrocyte 2' in KB-E0010.blend stack
    :param segment_range: A (start, end) tuple specifying which segment nodes to process
                          A range of (n,n+1) will create cut faces for all the points in the segment,
                          otherwise the cuts are created for each start segment node.
    '''

    dpath = '/var/remote/projects/epfl/data/KB-E0010/astrocyte2'
    skel_am_file = dpath + '/GeometrySurface.Smt.SptGraph_OK.am'
    skel_json_file = dpath+'/GeometrySurface.Smt.SptGraph_OK.annotations.json'
    obj_name = 'Astrocyte 2'

    reader = AmirameshReader()
    with open(skel_am_file, 'r') as f:
        skeleton = reader.parse(f)
    with open(skel_json_file, 'r') as f:
        data = json.load(f)

    debug_cut_faces(obj_name, skeleton, segment_range)


def test1():
    n_pos = Vector((17.5, 20.25, 11.5))
    n_pos = Vector((3.243425607681274, 0.7165750265121460, 24.52972412109375))
    n0_pos = Vector((3.143123149871826, 0.7475772500038147, 24.45019721984863))
    n_norm = n_pos - n0_pos

    # Create Cross Sections
    generate_node_cross_section_data('Astrocyte 2', n_pos, n_norm)


def test2():
    # import os, sys, imp, bmesh, functools
    # sys.path.append(os.path.abspath('/var/remote/projects/epfl/development/staging/skeletonizer/Source'))
    # import skeleton_annotate
    # skeleton_annotate.test2()
    # imp.reload(skeleton_annotate) # to reload changes
    dpath = '/var/remote/projects/epfl/data/KB-E0010/astrocyte2'
    am_path = dpath + '/GeometrySurface.Smt.SptGraph_OK.am'
    annotate_path = dpath+'/GeometrySurface.Smt.SptGraph_OK.annotations.json'
    segment_range = (0,1)
    generate_cross_sections('Astrocyte 2', am_path, annotate_path, segment_range)



def main():
    # blender -t 1 -b /var/remote/projects/epfl/data/KB-E0010/KB-E0010.blend -P skeleton_annotate.py -- "$CELLNAME" $AMPATH $START_SEG $SEG_SIZE

    argv = sys.argv[sys.argv.index("--") + 1:]

    # TODO: Add error handling
    cellname = argv[0]
    skel_path, skel_name, skel_am_file, skel_json_file = get_paths(argv[1])

    start_seg = max(0, int(argv[2]))
    end_seg = start_seg + max(1, int(argv[3]))
    segment_range = (start_seg, end_seg)

    generate_cross_sections(cellname, skel_am_file, skel_json_file, segment_range, skel_path)

if __name__ == '__main__':
    main()

