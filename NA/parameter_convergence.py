#!/usr/bin/env python
"""
Created on Feb 22, 2013

@author: behry
"""
import sys
sys.path.append('./')
from gpdcreport2mdl import *
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
import numpy as np
rcParams['figure.subplot.left'] = 0.1
rcParams['figure.subplot.right'] = 0.95
rcParams['figure.subplot.top'] = 0.92
rcParams['figure.subplot.bottom'] = 0.1
rcParams['figure.subplot.wspace'] = 0.25
rcParams['figure.subplot.hspace'] = 0.32

savedir = '/home/behry/uni/data/results_pgc'
repfile = os.path.join(savedir, 'reports/canada_48.000000_-105.000000.report')
os.system('/usr/local/Geopsy.org/bin/gpdcreport ' + repfile + ' >tmp.mdl')
models, data, nd, ne, nrow, nlayer = get_models('tmp.mdl')
param = np.zeros((nlayer, ne + 1))
misfit = np.zeros(ne + 1)
prmidx = 2
# todo: why does ne and data.size differ?
for iter in xrange(ne + 1):
    for layer in xrange(nlayer):
        param[layer, iter] = models[nd * iter + prmidx * nlayer + layer]
    misfit[iter] = data[iter]
plt.figure(figsize=(10, 10))
for i in xrange(9):
    H, xedges, yedges = np.histogram2d(misfit, param[i, :], bins=(40, 25),
                                       normed=True)
    extent = (xedges[0], xedges[-1], yedges[0], yedges[-1])
    plt.subplot(3, 3, i + 1)
    plt.imshow(np.log10(H.T), extent=extent, interpolation='nearest',
               cmap=cm.gray_r, aspect='auto')
    plt.title('Layer %d' % (i + 1))
    if i > 5:
        plt.xlabel('RMS misfit')
    if i % 3 == 0:
        plt.ylabel('S-velocity [km/s]')
plt.savefig(os.path.join(savedir, 'plots/s_velocity_convergence.png'))
plt.show()
