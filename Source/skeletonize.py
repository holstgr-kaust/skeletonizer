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

    for segm in skel.segments :
        if len(segm.points) < 2 :
            continue

        start, end = segm.points[0:2]
        # Given that the soma is empty, the first segment start point is
        # placed at 0, 0, 0. The position is corrected after adding the
        # first section. However, we cannot correct the diamater.
        section = soma.grow(end.x, end.y, end.z, end.diameters[0],
                            Section_Type.DENDRITE)
        section.move_point(0, Vector3f(start.x, start.y, start.z))

        for pt in segm.points[2:] :
            section.grow(pt.x, pt.y, pt.z, pt.diameters[0])

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
