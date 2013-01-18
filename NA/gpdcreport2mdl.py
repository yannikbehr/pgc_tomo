#!/usr/local/bin/python
"""
read models from *.report files into python variables
"""
from numpy import *
import os, os.path, sys

def get_models(filename):
    """
    read in dinver-model output
    """
    #ndmax = 450
    #nemax = 41000
    
    f = open(filename,'r')
    ### read header
    #f.readline()
    #nemax = int(f.readline().split()[4])
    nemax = int(f.readline().split()[1])
    curpos = f.tell()
    line = f.readline()
    nlayer = int(f.readline())
    nd = 0
    nrow = 0
    ncol = 0
    while True:
        line = f.readline()
        if line.find('model') != -1: break
        a = line.split()
        nd = nd + len(a)
        ncol = len(a)
        nrow += 1
    f.seek(curpos)
    models = zeros(nd*nemax)
    data = zeros(nemax)
    ne = -1

    while True:
        line = f.readline()
        if not line: break
        if line.find('model'):
            tmp = line.split()
            ind = int(tmp[3][0])
            ne += 1
            data[ne] = float(tmp[4].split('=')[1])
            f.readline()
            for _i in xrange(nrow):
                line = f.readline()
                a = line.split()
                for _j in xrange(ncol):
                    models[ne*nd+_i+_j*nrow] = float(a[_j])/1000.

    return (models,data,nd,ne,nrow,nlayer)


if __name__ == '__main__':
    try:
        mdlname = sys.argv[1]
    except:
        print "usage: %s report-file"%os.path.basename(sys.argv[0])
        sys.exit(1)
    os.system('gpdcreport '+mdlname+' >tmp.mdl')
    models, data, nd, ne = get_models('tmp.mdl')
    os.remove('tmp.mdl')
    print "number of models:     ",ne+1
    print "number of parameters: ",nd
    
