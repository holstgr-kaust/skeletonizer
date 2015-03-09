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


class MorphologyFileTestCase(unittest.TestCase):
    options = MorphologyCreateOptions()
    test_dir_path = os.path.abspath(os.path.split(__file__)[0])
    data_dir_path = os.path.join(test_dir_path,'data')

    def setUp(self):
        k_FORMAT = "%(message)s" # "%(asctime)-15s %(message)s"
        logging.basicConfig(format=k_FORMAT, level=self.options.verbosity_level)

        source_path = os.path.join(self.data_dir_path, 'test.SptGraph')

        self.options.set_pathname(source_path)
        self.options.skel_out_path = self.test_dir_path
        self.options.scaling_factor = 20.0
        self.options.set_filepaths()

    def tearDown(self):
        logging.shutdown()

    def test_create_mophology_file(self):
        try:
            os.remove(self.options.skel_out_file)
        except OSError:
            pass

        self.assertEqual(os.path.exists(self.options.skel_out_file), False, 
                        ('expected no %s file' % self.options.skel_out_file))

        self.options.validate()

        reader = AmirameshReader()

        with open(self.options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(self.options.skel_json_file, 'r') as f:
            data = json.load(f)

        self.options.set_annotation_data(data)

        morphology = create_morphology(skel, data['soma'], self.options)

        create_morphology_file(morphology, self.options)

        self.assertEqual(os.path.exists(self.options.skel_out_file), True, 
                        ('expecting %s output file' % self.options.skel_out_file))


        """
        # NOTE: display_morphology_file requires a path, not just a filename.
        rtneuron-app.py
        display_morphology_file('./test/test.SptGraph.h5')
        """

suite = unittest.TestLoader().loadTestsFromTestCase(MorphologyFileTestCase)
unittest.TextTestRunner(verbosity=2).run(suite)

