#!/usr/bin/env python

"""
This peice of code reads in an amira mesh to prase its content into a python object.
The obtained object will be processed using the BBPsdk to form a morphology that is
output in an hd5 file. The created hd5 file is saved in current directory as source.
"""
import os
from bbp_import_module import *
from amiramesh import *

def create_morphology(skel):
    """
	creates morphology from the skeleton obtained
	:param skel: skeleton data structure from amiramesh reader
	:return: BBPsdk morphology of the skeleton
	"""

    morphology = Morphology()
    soma = morphology.soma()

    # TODO: fix - calculation of nonroot / root is a hack that doesn't capture loops
    nodes = {}
    nonroot_idxs = [segm.end for segm in skel.segments]
    root_idxs = [nidx for nidx in range(0,len(skel.nodes)) if nidx not in nonroot_idxs]

    for segm in skel.segments:
        if segm.start in nonroot_idxs:
            continue

    	rnode_data = segm.points[0]
        rnode_pos = (rnode_data.x,rnode_data.y,rnode_data.z)
        rnode = soma.grow(rnode_pos[0], rnode_pos[1], rnode_pos[2], rnode_data.diameters[0], 
                          Section_Type.DENDRITE)
        rnode.move_point(0, Vector3f(rnode_pos[0], rnode_pos[1], rnode_pos[2]))
        nodes[rnode_pos] = rnode
        print('Root Node:'+segm.start)

    for segm in skel.segments :
        if len(segm.points) < 2 :
            continue

        node_data, section_start = segm.points[0:2]
        node_pos = (node_data.x,node_data.y,node_data.z)
	
        #print('Segment Start:'+segm.start+' End:'+segm.end)

        if node_pos in nodes:
            node = nodes[node_pos]
        else: # TODO: fix hack to handle newly found roots, because root / nonroot logic is incomplete
            node = soma.grow(node_pos[0], node_pos[1], node_pos[2], node_data.diameters[0], 
                             Section_Type.DENDRITE)
            node.move_point(0, Vector3f(node_pos[0], node_pos[1], node_pos[2]))
            nodes[node_pos] = node
            print('New Root Node:'+segm.start)

        section = node.grow(section_start.x, section_start.y, section_start.z, section_start.diameters[0],
                            Section_Type.DENDRITE)

        # TODO: Sections should end at a node (unless there are loops)
        for pt in segm.points[2:] :
            sec_pos = (pt.x, pt.y, pt.z)
            if (sec_pos in nodes):
                section.grow(pt.x, pt.y, pt.z, pt.diameters[0])
                # TODO: if (sec_pos in nodes) section.grow(nodes[sec_pos]) else
                #section.grow(nodes[sec_pos])
            else:
                section.grow(pt.x, pt.y, pt.z, pt.diameters[0])

        node_data = segm.points[-1:][0]
        node_pos = (node_data.x,node_data.y,node_data.z)
        if node_pos not in nodes:
            nodes[node_pos] = section.grow(node_pos[0], node_pos[1], node_pos[2], node_data.diameters[0], 
                                Section_Type.DENDRITE)
            #print('New Node:'+segm.end)

    return morphology


def create_morphology_file(label, skel):
    """
    takes the morphology object after calling create_morphology function and creates and hdf5 file for it
    :param label: name of the amiramesh file in which hdf5 file is named after
    """

    try:
        os.remove(label + '.h5')
    except OSError:
        pass

    morphology = create_morphology(skel)
    morphology.label(label)

    try:
        directory = os.path.abspath(os.path.curdir)
        writer = Morphology_Writer()
        writer.open(directory)
        writer.write(morphology, Morphology_Repair_Stage.RAW_MORPHOLOGY)
    except OSError:
        pass


if __name__ == '__main__':
    if(len(sys.argv)==1):
        print("amira mesh file is required")
    else:
        label = sys.argv[1]
        print ("HDF5 Skeletonizer")
        print label
        reader = AmirameshReader()

        with open(str(label), 'r') as f:
            skel = reader.parse(f)

        create_morphology_file(label, skel)

    """
    rtneuron-app.py
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.juan.h5')
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.h5')
    """
