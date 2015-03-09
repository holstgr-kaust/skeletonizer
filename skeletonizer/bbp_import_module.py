#!/usr/bin/python

import sys
import os

try:
    import bbp
except ImportError:
    sys.path = [  os.path.expanduser('~')+'/Development/RTNeuron/Build/BBPSDK/lib'
                , '/var/remote/projects/epfl/development/staging/RTNeuron/Build/BBPSDK/lib'
                , '/var/remote/projects/epfl/development/production/RTNeuron/Build/BBPSDK/lib'
               ] + sys.path

from bbp import *
