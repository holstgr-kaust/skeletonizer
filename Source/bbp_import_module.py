#!/usr/bin/python

import sys
import os

sys.path = ['/var/remote/projects/epfl/development/staging/RTNeuron/Build/BBPSDK/lib',
            os.path.expanduser('~')+'/Development/RTNeuron/Build/BBPSDK/lib'] + sys.path

from bbp import *
