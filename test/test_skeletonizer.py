"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
    Skeletonizer testing module 

    To run tests:
    python -m unittest test.test_skeletonizer
"""

import unittest
import os
import re
import sys
import math
import getopt
import copy
import csv
import json
import logging
import operator
from collections import defaultdict
import subprocess

try:
    import skeletonizer
except ImportError:
    sys.path.append(os.path.abspath(os.path.dirname(os.path.abspath(os.path.split(__file__)[0]))))

from skeletonizer.bbp_import_module import *
from skeletonizer.amiramesh import *
from skeletonizer.maths import *
from skeletonizer.graphs import *
from skeletonizer.morphology import *


class MorphologyFileTestCase(unittest.TestCase):
    test_dir_path = os.path.abspath(os.path.split(__file__)[0])
    data_dir_path = os.path.join(test_dir_path, 'data')
    bin_dir_path = os.path.join(os.path.abspath(os.path.dirname(test_dir_path)), 'bin')

    def setUp(self):
        k_FORMAT = "%(message)s" # "%(asctime)-15s %(message)s"
        logging.basicConfig(format=k_FORMAT, level=logging.WARNING)

    def tearDown(self):
        logging.shutdown()

    def test_create_mophology_file(self):
        source_path = os.path.join(self.data_dir_path, 'test.SptGraph')

        options = MorphologyCreateOptions()
        options.set_pathname(source_path)
        options.skel_out_path = self.test_dir_path
        options.scaling_factor = 20.0
        options.set_filepaths()

        try:
            os.remove(options.skel_out_file)
        except OSError:
            pass

        self.assertEqual(os.path.exists(options.skel_out_file), False, 
                        ('expected no %s file' % options.skel_out_file))

        options.validate()

        reader = AmirameshReader()

        with open(options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(options.skel_json_file, 'r') as f:
            annotation_data = json.load(f)

        options.set_annotation_data(annotation_data)


        with open(options.skel_csv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', quotechar='|')
            xsection_data = {}
            for r in reader:
                area = float(r['area'])
                perimeter = float(r['perimeter'])
                diameter = math.sqrt(area) / math.pi
                xsection_data[(int(r['segment_idx']), int(r['pnt_idx']))] = \
                            {'area':area, 'perimeter':perimeter, \
                             'estimated_diameter':float(r['estimated_diameter']), \
                             'estimated_area':float(r['estimated_area']), \
                             'estimated_perimeter':float(r['estimated_perimeter']), \
                             'blender_position':r['blender_position'], \
                             'blender_normal':r['blender_normal'], \
                             'diameter':diameter}

        skel.update_diameters(xsection_data)
        options.set_xsection_data(xsection_data)


        morphology = create_morphology(skel, annotation_data['soma'], options)

        create_morphology_file(morphology, options)

        self.assertEqual(os.path.exists(options.skel_out_file), True, 
                        ('expecting %s output file' % options.skel_out_file))

        # Open created morphology file in RTNeuron to view
        # rtneuron-app.py doesn't support this yet
        with open(os.devnull, 'w') as nof:
            if not subprocess.call(['which','rtneuron-app.py'], stdout=nof, stderr=nof):
                #subprocess.call(['rtneuron-app.py', '--script', 'display_morphology_file(%s)' % options.skel_out_file]) 
                pass

    def test_create_crosssection_file(self):
        source_path = os.path.join(self.data_dir_path, 'test.SptGraph')
        blend_file = os.path.join(self.data_dir_path, 'test.blend')
        script_file = os.path.join(self.bin_dir_path, 'skeleton_annotate_csv.py')

        self.assertEqual(os.path.exists(script_file), True, 
                        ('expecting %s script file' % script_file))
        self.assertEqual(os.path.exists(blend_file), True, 
                        ('expecting %s input file' % blend_file))

        options = MorphologyCreateOptions()
        options.set_pathname(source_path)
        options.skel_out_path = self.test_dir_path
        options.force_overwrite = True # we're not writing the *.h5 morphology file, so ignore when validating
        options.set_filepaths()

        options.validate()

        reader = AmirameshReader()

        with open(options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(options.skel_json_file, 'r') as f:
            data = json.load(f)

        options.set_annotation_data(data)

        r = (0, len(skel.segments))

        # TODO: fix ugly duplicate name hack -- we should specify the name of the output file
        obj_name = 'shape'
        csvfilename = obj_name+'-cross_section_data-range-'+str(r[0])+'-'+str(r[1]-1)+'-of-'+str(len(skel.segments))
        csv_out_file = os.path.join(self.data_dir_path, csvfilename + '.csv')
        
        try:
            os.remove(csv_out_file)
        except OSError:
            pass

        self.assertEqual(os.path.exists(csv_out_file), False, 
                        ('expected no %s file' % csv_out_file))

        # invoke blender with script
        #

        with open(os.devnull, 'w') as nof:
            self.assertEqual(subprocess.call(['which','blender'], stdout=nof, stderr=nof), 0, 
                        'expected to find blender executable file')

        subprocess.call(['blender', '-b', blend_file, '-P', script_file, '--', obj_name, options.skel_am_file, str(r[0]), str(r[1])]) 

        self.assertEqual(os.path.exists(csv_out_file), True, 
                        ('expecting %s output file' % csv_out_file))

        # TODO: Scan stdout from subprocess.call to find errors or issues (e.g., "No cross-section data for node:")


suite = unittest.TestLoader().loadTestsFromTestCase(MorphologyFileTestCase)
unittest.TextTestRunner(verbosity=2).run(suite)

