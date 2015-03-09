"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
    Skeletonizer testing module 
"""

import unittest
import os
import re
import sys
import math
import getopt
import copy
import json
import logging
import operator
from collections import defaultdict

try:
    import skeletonizer
except ImportError:
    sys.path.append(os.path.abspath(os.path.dirname(os.path.abspath(os.path.split(__file__)[0]))))

from skeletonizer.bbp_import_module import *
from skeletonizer.amiramesh import *
from skeletonizer.maths import *
from skeletonizer.graphs import *
from skeletonizer.morphology import *


class MorphologyTestCase(unittest.TestCase):
    options = MorphologyCreateOptions()
    test_dir_path = os.path.abspath(os.path.split(__file__)[0])
    data_dir_path = os.path.join(options.skel_path,'data')
    source_path = os.path.join(data_dir_path, 'test.SptGraph')

    def setUp(self):
        k_FORMAT = "%(message)s" # "%(asctime)-15s %(message)s"
        logging.basicConfig(format=k_FORMAT, level=options.verbosity_level)

        self.options.set_pathname(source_path)
        self.options.skel_out_path = self.test_dir_path
        self.options.scaling_factor = 20.0
        self.options.set_filepaths()

    def tearDown(self):
        logging.shutdown()

    def test_create_mophology_file(self):
        output_file = 'test.SptGraph.h5'

        print("options: %s " % (self.options.skel_path,self.options.skel_name,self.options.skel_out_path,
                                self.options.skel_am_file,self.options.skel_json_file,self.options.skel_out_file))
        '''
        try:
            os.remove(self.options.skel_out_file)
        except OSError:
            pass

        self.assertEqual(os.path.exists(self.options.skel_out_file), False, 
                        ('expected no %s file' % self.options.skel_out_file))

        self.options.validate()

        reader = AmirameshReader()

        with open(options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(options.skel_json_file, 'r') as f:
            data = json.load(f)

        self.options.set_annotation_data(data)

        morphology = create_morphology(skel, data['soma'], self.options)

        create_morphology_file(morphology, options)
        '''
        self.assertEqual(os.path.exists(self.options.skel_out_file), True, 
                        ('expecting %s output file' % self.options.skel_out_file))


    """
    # NOTE: display_morphology_file requires a path, not just a filename.
    rtneuron-app.py
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.juan.h5')
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.h5')
    """
