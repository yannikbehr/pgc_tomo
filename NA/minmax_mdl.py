#!/usr/bin/env python
"""
Plot a minimum and maximum model based on the
dinver parameterisation file.
"""

from gpdcparam2model import ProFile
import sys
sys.path.append('../3D')
from gmtpy import GMT
import matplotlib.pyplot as plt
import numpy as np
def minmax_mdl(fname):
    prf = ProFile(fname)
    minmdl = []
    maxmdl = []
    ymin = 0
    ymax = 0
    for i in xrange(prf.nlayers):
        print prf.vsmin[i], prf.vsmax[i]
        pmin1 = [prf.vsmin[i], ymax]
        pmin2 = [prf.vsmin[i], ymax - prf.thmax[i]]
        pmax1 = [prf.vsmax[i], ymin]
        pmax2 = [prf.vsmax[i], ymin - prf.thmin[i]]
        ymax -= prf.thmax[i]
        ymin -= prf.thmin[i]
        minmdl.append(pmin1)
        minmdl.append(pmin2)
        maxmdl.append(pmax1)
        maxmdl.append(pmax2)
    # add the half space
    yhalf = ymax - 50.  # just a random number
    pminhalf = [prf.vsmin[i], yhalf]
    minmdl.append(pminhalf)
    pmaxhalf = [prf.vsmax[i], yhalf]
    maxmdl.append(pmaxhalf)
    return np.array(minmdl), np.array(maxmdl)


if __name__ == '__main__':
    fname = '/home/behry/uni/data/results_pgc/param/Lith50_Canada_48.000000_-105.000000.param'
    mdl_min, mdl_max = minmax_mdl(fname)
    print mdl_min
    print mdl_max
    plt.plot(mdl_min[:, 0], mdl_min[:, 1], 'b-')
    plt.plot(mdl_max[:, 0], mdl_max[:, 1], 'r-')
    plt.show()
