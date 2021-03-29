from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy

import landlab
from landlab import RasterModelGrid, FIXED_VALUE_BOUNDARY
from landlab.components import (FlowRouter,
                                FastscapeEroder,
                                LinearDiffuser)
from landlab.components.flexure import Flexure
#from landlab.components.gflex import gFlex
from flexure import gFlex
from landlab.plot.imshow import imshow_grid
from landlab.plot.drainage_plot import drainage_plot
from landlab.io.netcdf import write_netcdf
from landlab.io.netcdf import write_raster_netcdf
from landlab.io.netcdf import read_netcdf

sys.path.append(os.environ['PARAMPARSER_PATH'])
from paramparser import ParamParser
from paramparser import MissingKeyError

from datasaver import DataSaver
