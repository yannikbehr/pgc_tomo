#!/usr/bin/env mypython
"""
Plot Moho as a 2D surface 
"""

from enthought.mayavi import mlab
mlab.options.offscreen = True
from enthought.tvtk.api import tvtk
from enthought.mayavi.sources.api import VTKDataSource
from scipy.io import savemat,loadmat
from collections import defaultdict
from gmtpy import GMT
import pdb
import sys
import os
import numpy as np
sys.path.append('../NA')
from plot_1D import get_best_model
from plot_3D import get_1Dresults, get_tvtk_grid
from gpdcreport2mdl import *
from dplot import dplot
import tempfile
import delaz
import cStringIO
