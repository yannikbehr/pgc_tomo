#!/usr/bin/env python
"""
Generate synthetic dispersion curve for a given 1D earth model.
"""

import os
import sys
import tempfile
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt

def get_disp(model, fmin, fmax, gv=False, nfreq=100, nmode=1, gpdcbin='/usr/local/Geopsy.org/bin/gpdc'):
    """
    Calculate synthetic dispersion curve using geopsy's gpdc. On success it returns
    two 2D Numpy arrays. The first one for Rayleigh waves the second for Love waves.
    model: string containing the 1D earth model
    fmin: minimum frequency of the dispersion curves
    fmax: maximum frequency of the dispersion curves
    gv: if True group velocity is returned, otherwise phase velocity
    nfreq: number of frequency samples
    nmode: number of modes
    gpdcbin: path to local gpdc executable
    return: The dimensions of the returned arrays are 2*nfreqx2*nmode.
            The dispersion curves are returned in slowness vs frequency
    """
    mfile = tempfile.mktemp()
    dfile = tempfile.mktemp()
    if os.path.isfile(mfile):
        os.remove(mfile)
    if os.path.isfile(dfile):
        os.remove(dfile)
    open(mfile, 'w').writelines(model)
    if gv:
        cmd = "%(gpdcbin)s -R %(nmode)d -L %(nmode)d -min %(fmin)f -max %(fmax)f -group -n %(nfreq)d < %(mfile)s > %(dfile)s" % vars()
    else:
        cmd = "%(gpdcbin)s -R %(nmode)d -L %(nmode)d -min %(fmin)f -max %(fmax)f -n %(nfreq)d < %(mfile)s > %(dfile)s" % vars()
    os.system(cmd)
    rw = np.zeros((nfreq * 2, 2 * nmode))
    lw = np.zeros((nfreq * 2, 2 * nmode))
    f = open(dfile)
    while True:
        line = f.readline()
        if not line: break
        if line.find('Rayleigh') != -1:
            rcnt = 0
            ### jump the next two lines
            line = f.readline()
            line = f.readline()
            mode = 0
            ### now read in till the next comment
            while True:
                line = f.readline()
                if line.find('Mode') != -1:
                    mode += 2
                    rcnt = 0
                    line = f.readline()
                if line.find('Love') != -1: break
                if not line: break
                fq, s = map(float, line.split())
                rw[rcnt, mode] = fq
                rw[rcnt, mode + 1] = s
                rcnt += 1
        if line.find('Love') != -1:
            lcnt = 0
            ### jump the next two lines
            line = f.readline()
            line = f.readline()
            mode = 0
            ### now read in till the next comment
            while True:
                line = f.readline()
                if line.find('Mode') != -1:
                    mode += 2
                    lcnt = 0
                    line = f.readline()
                if not line: break
                fq, s = map(float, line.split())
                lw[lcnt, mode] = fq
                lw[lcnt, mode + 1] = s
                lcnt += 1

    os.remove(mfile)
    #os.remove(dfile)
    return rw, lw

if __name__ == '__main__':
    plotdisp = True
    gpdcbin = '/usr/local/Geopsy.org/bin/gpdc'
    ### first create synthetic dispersion curve
    ### crust 2.0 model for CNIPSE + upper mantle from PREM
    #number of layers including halfspace
    #1st layer: Thickness [m] Vp [m/s] Vs [m/s] Density [kg/m3]
    #2nd layer: ....
    #halfspace: Thickness=0.
    model = """
8
1000.0  2500.0    1200.0    2100.0
1000.0  4000.0    2100.0    2400.0
16000.0 6000.0    3500.0    2700.0
8000.0  6600.0    3700.0    2900.0
9000.0  7200.0    4000.0    3050.0
20000.0 8079.0    4465.0    3377.0
20000.0 8067.0    4457.0    3375.0
0.0     8067.0    4457.0    3375.0
"""

    rayc, lovc = get_disp(model, 0.02, 1.0, nmode=2, gpdcbin=gpdcbin)
    rayu, lovu = get_disp(model, 0.02, 1.0, gv=True, gpdcbin=gpdcbin)
    indc = np.where(lovc[:, 0] > 0.)
    indu = np.where(lovu[:, 0] > 0.)
    if plotdisp:
        plt.plot(1. / lovc[indc, 0][0], 1. / lovc[indc, 1][0] / 1000., label='Love phase')
        plt.plot(1. / rayc[indc, 0][0], 1. / rayc[indc, 1][0] / 1000., label='Rayleigh phase')
        plt.plot(1. / lovu[indu, 0][0], 1. / lovu[indu, 1][0] / 1000., label='Love group')
        plt.plot(1. / rayu[indu, 0][0], 1. / rayu[indu, 1][0] / 1000., label='Rayleigh group')
        plt.xlabel('Period [s]')
        plt.ylabel('Velocity [km/s]')
        plt.legend(loc='lower right')
        plt.show()




