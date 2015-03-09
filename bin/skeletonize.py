#!/usr/bin/env python

"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
This program creates a BBPSDK HDF5 morphology file from a skeletonization representation in an Amiramesh text file.
"""

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


if __name__ == '__main__':
    options = MorphologyCreateOptions()
    k_FORMAT = "%(message)s" # "%(asctime)-15s %(message)s"
    logging.basicConfig(format=k_FORMAT, level=options.verbosity_level)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hfas:o:v:t:x:",["skeleton=","output_dir=","verbose=","threshold=","scale="])
    except getopt.GetoptError:
        print 'skeletonize.py -h'
        sys.exit(2)
    else:
        for opt, arg in opts:
            if opt == '-h':
                print 'Skeletonize converts an Amiramesh skeleton graph, plus annotations, into a BBPSDK cell morphology.'
                print '\nUsage:'
                print ' skeletonize.py <skeleton>'
                print ' skeletonize.py [-v <level>] [-a] [-t <threshold>] [-x <scale>] -s <skeleton> [-f] [-o <output_dir>]'
                print '\t -a \t\t Allow cycles in skeleton graph (default False)'
                print '\t -f \t\t Force overwrite of output files'
                print '\t -o <dirname>\t Output directory'
                print '\t -s <filename>\t Input skeleton filename'
                print '\t -t <threshold>\t Set minimum segment arc length (default 0)'
                print '\t -v <level>\t Set verbosity level: %i-%i' % (logging.NOTSET, logging.FATAL)
                print '\t -x <scale>\t Set skeleton scaling factor to resize output skeleton'
                print '\nExample:'
                print '\t # creates /<path>/cell.Smt.SptGraph.h5 from /<path>/cell.Smt.SptGraph'
                print '\t skeletonize.py -s cell.Smt.SptGraph'
                print '\nNotes:'
                print '\t For input source <filename>, expected input files are:'
                print '\t\t <filename>.am # Amiramesh text file of skeleton graph'
                print '\t\t <filename>.annotations.json # JSON file with {"soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r}}'
                print '\t\t\t Measurements such as "centre" and "radius" are in the coordinate system and units of the input source.'
                print '\t Output file(s) are:'
                print '\t\t <filename>.h5 # BBPSDK HDF5 format'
                print '\t Verbosity levels(s) are:'
                print '\t\t all=0, debug=10, INFO=20, warning=30, error=40'
                print '\t\t INFO is the default logging level'
                print '\t\t Debug logging levels include visual debugging artifacts added to the morphology file.'
                print '\t\t Visual debugging includes:'
                print '\t\t\t Soma star: representation of soma size and location.'
                print '\t\t\t Coordinate axis: X, Y, Z are represented as three bars with end-fingers (0=X,1=Y,2=Z).'
                print '\t\t All logging level includes additional visual debugging artifacts:'
                print '\t\t\t Soma dendrites: visual representation of original source soma skeleton.'
                print '\t Scale specifies the final scaling factor applied to the output files.'
                print '\t Threshold specifies the minimum segment section length in original unscaled units.'
                print '\t Display in rtneuron-app.py using: display_morphology_file(\'/<path>/<filename>.h5\')'
                print '\t\t NOTE: \'display_morphology_file\' may require a relative or absolute path, not just a filename, to display morphology.'
                sys.exit()
            elif opt == '-a':
                options.allow_cycles = True
                logging.info("Allow Cycles set to: %s", options.allow_cycles)
            elif opt == '-f':
                options.force_overwrite = True
            elif opt in ("-o", "--output_dir"):
                options.skel_out_path = arg
                if (not os.path.isdir(options.skel_out_path)):
                    logging.error('ERROR - Output directory must be directory:%s', options.skel_out_path)
                    sys.exit(4)
            elif opt in ("-s", "--skeleton"):
                options.set_pathname(arg)
            elif opt in ('-t', "--threshold"):
                options.force_segment_threshold = True
                options.threshold_segment_length = float(arg)
                logging.info("Segment length threshold set to: %f", options.threshold_segment_length)
            elif opt in ('-v', "--verbose"):
                options.verbosity_level = int(arg)
                logging.getLogger().setLevel(options.verbosity_level)
                logging.info("Verbosity set to: %i", options.verbosity_level)
            elif opt in ('-x', "--scale"):
                options.scaling_factor = float(arg)
                logging.info("Morphology scaling factor set to: %f", options.scaling_factor)

        if not opts:
            if len(sys.argv) != 2:
                logging.error('Expected skeleton source. Try: skeletonize.py -h')
                sys.exit(2)
            options.set_pathname(sys.argv[1])

        options.set_filepaths()

        options.validate()

        logging.info('HDF5 Skeletonizer')
        logging.info('\t Source graph: %s', options.skel_am_file)
        logging.info('\t Source annotations: %s', options.skel_json_file)
        if options.force_overwrite:
            logging.info('\nFORCING OVERWRITE of output file: %s\n', options.skel_out_file)

        reader = AmirameshReader()

        with open(options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(options.skel_json_file, 'r') as f:
            data = json.load(f)

        options.set_annotation_data(data)

        morphology = create_morphology(skel, data['soma'], options)

        create_morphology_file(morphology, options)

        logging.info('Wrote out file: %s', options.skel_out_file)

    finally:
        logging.shutdown()

    """
    # NOTE: display_morphology_file requires a path, not just a filename.
    rtneuron-app.py
    display_morphology_file('./GeometrySurface.Smt.SptGraph.h5')
    """
