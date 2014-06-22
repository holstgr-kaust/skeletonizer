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

    for segm in skel.segments:

        print("Segment start:" + repr(segm.start) + "end:" + repr(segm.end))
        prev_node = soma
        prev_diameter = None
        for pt in segm.points:
            prev_diameter = prev_diameter if prev_diameter else pt.diameters[0]

            next_node = prev_node.grow(pt.x, pt.y, pt.z, pt.diameters[0], Section_Type.DENDRITE, UNDEFINED_COUNT())

            prev_diameter = pt.diameters[0]
            prev_node = next_node

            ## Print Out
            ptstr = "\t<" + repr(pt.x) + "," + repr(pt.y) + "," + repr(pt.z) + ">["
            ##diametere loop
            for diam in pt.diameters:
                ptstr += repr(diam) + ","
            ptstr += "]"
            print(ptstr)
        print ""

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
