# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# ***** END GPL LICENCE BLOCK *****

bl_info = {
    "name": "Cross Section",
    "author": "Yorik van Havre, Alejandro Sierra, Howard Trickey, Campbell Barton, Ejnar Rasmussen",
    "description": "Creates cross section(s) of the selected object(s) using the active object as cut plane",
    "version": (0, 1, 8),
    "blender": (2, 6, 5),
    "category": "Object",
    "location": "Toolshelf > Cross Section",
    "warning": '',
    "wiki_url": "",
    "tracker_url": "",
    }


"""Script made by Yorik for 2.4x based upon a code by Cambo found on blenderArtists, 
and ported to 2.5x by Ejnaren. 

This scripts creates cross-sections of selected objects, at intersection with 
active plane object. 

Installation:

Via the addons menu, no particular action needed

Usage:

Select objects you want to cut, then select (make active) 
the cutting plane, then run the script. The resulting section parts will be 
filled and nicely placed inside a group for easy selection.

Options:

You can turn the fill option on or off, if enabled, closed edge loops are turned into faces.

Limitations:

- Only Mesh and Surface (by approximation) objects will be cut
- The cutting object must be a plane (it must have only one face)
- The cutting plane shouldn't have any parents, if it does, the rotation from 
  the parents will not affect the section's position and rotation"""

import bpy, threading, time
from mathutils import *
from math import *


def centerMass(me):
    """Finds the Center of Mass on a mesh ( no Vertex Weights taken into account ).   
    me: Blender Mesh - the mesh whose center you want to find
    Returns: Vector -  The center of mass of the given geometry
    NOTE: The Center of Mass differs from the Bounding Box Center and is quite 
    usefull for 2D Shapes."""
    
    sum = Vector()
    for v in me.vertices:
        sum += v.co
    center = (sum)/len(me.vertices)
    return center


def invRotation(mx):
    """Finds the matrix that does the invert rotation of the given matrix but takes
    scale out of the way to keep things sane
    mx: Matrix - the matrix whose rotation you want to invert
    Returns: Matrix - a matrix of the inverse rotation"""
    
    r = mx.to_4x4()
    s = mx.to_scale()
    smx = Matrix()
    for i in range(0,3):
        smx[i][i] = s[i]
    rmx = r.inverted() * smx
    return rmx


def dupTest(object):
    """Checks objects for duplicates enabled (any type)
    object: Blender Object.
    Returns: Boolean - True if object has any kind of duplicates enabled."""
    
    if (object.is_duplicator):
       return True
    else:
       return False


def getObjectsAndDuplis(oblist,MATRICES=False,HACK=False):
    """Return a list of real objects and duplicates and optionally their matrices
    oblist: List of Blender Objects
    MATRICES: Boolean - Check to also get the objects matrices default=False
    HACK: Boolean - See note default=False
    Returns: List of objects or
             List of tuples of the form:(ob,matrix) if MATRICES is set to True  
    NOTE: There is an ugly hack here that excludes all objects whose name 
    starts with "dpl_" to exclude objects that are parented to a duplicating 
    object, User must name objects properly if hack is used."""
    
    result = []
        
    for ob in oblist:
        if dupTest(ob):
            dup_obs=ob.children
            if len(dup_obs):
                for dup_ob in dup_obs: #, dup_mx
                    if MATRICES:                        
                        result.append((dup_ob,dup_ob.matrix_world))
                    else:
                        result.append(dup_ob)
        else:
            if HACK:
                if ob.name[0:4] != "dpl_":
                    if MATRICES:
                        mx = ob.matrix_world
                        result.append((ob,mx))
                    else:
                        result.append(ob)
            else:
                if MATRICES:
                    mx = ob.matrix_world
                    result.append((ob,mx))
                else:
                    result.append(ob)
            
    return result


def section(cut_me,mx,pp,pno,FILL=True):
    """Finds the section mesh between a mesh and a plane 
    cut_me: Blender Mesh - the mesh to be cut
    mx: Matrix - The matrix of object of the mesh for correct coordinates
    pp: Vector - A point on the plane
    pno: Vector - The cutting plane's normal
    FILL: Boolean - Check if you want to fill the resulting mesh, default=True
    Returns: Mesh - the resulting mesh of the section if any or
             Boolean - False if no section exists""" 

    verts = []
    ed_xsect = {}

    for ed in cut_me.edges:

        # getting a vector from each edge vertices to a point on the plane  
        # first apply transformation matrix so we get the real section
        
        vert1 = ed.vertices[0]
        v1 = cut_me.vertices[vert1].co * mx.transposed() 
        co1 = v1 - pp
        vert2 = ed.vertices[1]
        v2 = cut_me.vertices[vert2].co * mx.transposed()
        co2 = v2 - pp

        # projecting them on the normal vector
        proj1 = co1.project(pno).length
        proj2 = co2.project(pno).length

        if (proj1 != 0):
            rad1 = co1.angle(pno)
            angle1 = round((180 * rad1 / pi),4)

        else: angle1 = 0
        if (proj2 != 0):
            rad2 = co2.angle(pno)
            angle2 = round((180 * rad2 / pi),4)

        else: angle2 = 0
        
        #Check to see if edge intersects. Also check if edge is coplanar to the
        #cutting plane (proj1=proj2=0)

        if ((proj1 == 0) or (proj2 == 0) or  \
            (angle1 > 90) != (angle2 > 90)) and  \
            (proj1+proj2 > 0):

            #edge intersects.
            
            proj1 /= proj1+proj2
            co = ((v2-v1)*proj1)+v1
            verts.append(co)

            #store a mapping between the new vertices and the mesh's edges
            ed_xsect[ed.key] = len(ed_xsect)

    edges = []
    for f in cut_me.polygons:
        # get the edges that the intersecting points form
        # to explain this better:
        # If a face has an edge that is proven to be crossed then use the
        # mapping we created earlier to connect the edges properly
        ps = [ ed_xsect[key] for key in f.edge_keys if key in ed_xsect]

        if len(ps) == 2:
            edges.append(tuple(ps))

    if edges:
    
        x_me = bpy.data.meshes.new('Section')
        x_me.from_pydata(verts,edges,[])
        
        #create a temp object and link it to the current scene to be able to 
        #apply rem Doubles and fill 
        tmp_ob = bpy.data.objects.new('Mesh', x_me)

        sce = bpy.context.scene
        sce.objects.link(tmp_ob)
        
        # do a remove doubles to cleanup the mesh, this is needed when there
        # is one or more edges coplanar to the plane.
        bpy.context.scene.objects.active = tmp_ob

        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.select_mode(type="EDGE", action="ENABLE")
        bpy.ops.mesh.select_all(action="SELECT")

        # remove doubles:
        bpy.ops.mesh.remove_doubles()

        if FILL:
            bpy.ops.mesh.edge_face_add()
        # recalculate outside normals:
        bpy.ops.mesh.normals_make_consistent(inside=False)

        bpy.ops.object.mode_set(mode='OBJECT')
        
        #Cleanup
        sce.objects.unlink(tmp_ob)
        del tmp_ob
        
        return x_me
    else:
        return False

# operator definition

class OBJECT_OT_cross_section(bpy.types.Operator):
    """Creates cross section(s) of the selected object(s), using the active object as cut plane"""
    bl_idname = "object.cross_section"
    bl_label = "Cross Section"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Creates cross section(s) of the selected object(s), using the active object as cut plane"
    
    @classmethod
    def poll(cls, context):
        if context.active_object:
            if len(context.selected_objects) > 1:
                return True
        return False

    def execute(self, context):

        sce = context.scene
        ob_act = context.active_object   
        sel_group = context.selected_objects
    
        if not ob_act or ob_act.type != 'MESH':
            self.report({'WARNING'}, "No meshes or active object selected")
            return
    
        if len(sel_group)>=2:
    
            #Window.WaitCursor(1) 2.6 equivalent?
            t = time.time()
    
            #we need to be in edge or vertex mode otherwise fill() won't work
            #global MODE
            #MODE = Mesh.Mode()
            #Mesh.Mode(Mesh.SelectModes.EDGE) 
            cp = ob_act.data
        
            #Get section Plane's transformation matrix and euler rotation       
            p_mx = ob_act.matrix_world
            p_rot = ob_act.rotation_euler
            
            #Get the plane normal vector and a point on the plane
            pno = cp.polygons[0].normal * p_mx.to_3x3().transposed()
    
            pp_rot = cp.vertices[0].co * p_mx.transposed()
            pp = cp.vertices[0].co * p_mx.transposed()
    
            #filter selection to get back any duplis to cut them as well
    
            oblist = getObjectsAndDuplis(sel_group,MATRICES=True,HACK=False)
            
            #deselect all selected objects so we can select new ones and put 
            #them into the group that will contain all new objects
            for o in sel_group:
                o.select=False
            
            #create a list to hold all objects that we'll create
            parts = []
            for o,mx in oblist:
                typ = o.type
                if o != ob_act and \
                (typ =="MESH" or typ=="SURFACE" or typ=="CURVE"):
                    
                    #Use BpyMesh to get mesh data so that modifiers are applied
                    #and to be able to cut surfaces and curves (curves behave 
                    #strange though and don't seem to be very usefull
    
                    cut_me = o.to_mesh(scene=context.scene, apply_modifiers=True, settings='PREVIEW')
                    
                    #Run the main function
                    x_me = section(cut_me,mx,pp,pno,context.scene.cross_section_fill)
                    
                    #if there's no intersection just skip the object creation 
                    if x_me:
                        
                        #x_me = bpy.data.meshes.new('partofsection')
                        
                        #part_me = bpy.data.meshes.new('partofsection')
                        part_ob = bpy.data.objects.new("Partofsection", x_me)       
    
                        # Reset section part's center and orientation
                        # Since it's a 2d object it's best to use the objects 
                        # center of mass instead of bounding box center
                        loc = part_ob.location
                        center = centerMass(x_me) + loc
                
                        lmx =  Matrix.Translation(loc-center)
                        
                        x_me.transform(lmx)
    
                        part_ob.location = center
                        
                        parts.append(part_ob)
                        sce.objects.link(part_ob)
                        
                        # select the parts so that the user can do sth with them
                        # immidiately after the script is done
                        
                        part_ob.select=True
                        context.scene.objects.active = part_ob
                        
                    else:
                        self.report({'WARNING'}, "an object does not intersect the plane, continuing happily")
            #Put the parts of the section into a new group so that it's easy to
            # select them all
            sect_grp = bpy.data.groups.new('section')
            for part in parts:
                sect_grp.objects.link(part)
            
        else: 
            self.report({'WARNING'}, "the selection is empty, no object to cut!")
        
        print ('CrossSection finished in %.2f seconds' % (time.time()-t) )
        
        #Window.WaitCursor(0) #2.6 equivalent? couldn't find any...
        
        ob_act.select = False   
        context.scene.update()        
            
        return {'FINISHED'}
    
# a panel containing 2 buttons
class VIEW3D_PT_tools_cross_section(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_label = "Cross Section"
    #bl_context = "object" # for some reason this prevents the tool panel from appearing

    def draw(self, context):
        row = self.layout.row(align=True)
        row.alignment = 'LEFT'
        row.prop(context.scene, "cross_section_fill")
        row.operator("object.cross_section", text="Create cross section")
        
# Registers the operator, the toolshelf panel and the fill property
def register():
    # this registers all the panels and operators
    bpy.utils.register_module(__name__)
    # this stores the fill property in the current scene
    bpy.types.Scene.cross_section_fill = bpy.props.BoolProperty(
            name="Fill",
            description="Fill closed contours with faces",
            default=True)

# Removes the operator, the toolshelf and the fill property
def unregister():
        bpy.utils.unregister_module(__name__)
        del bpy.types.Scene.cross_section_fill

# This lets you import the script without running it
if __name__ == "__main__":
    register()
