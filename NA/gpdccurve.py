#!/usr/local/bin/python

"""
plot best-fitting dispersion curve
"""

import sys, os, os.path
from subprocess import *
from pylab import *
from optparse import OptionParser

def gpdccurve(mdlfile,wtype='R',ptype='p',mode=0):
    p1 = Popen(["gpdcreport","-pm",mdlfile],stdout=PIPE)
    mf = 1e99
    while True:
        line = p1.stdout.readline()
        if not line: break
        if line.find('#') != -1: continue
        a = line.split()
        _m = float(a[int(a[1])+2])
        if _m < mf:
            mf = _m
            idx = a[0]
#    if opts.verbose:
#        print "index of best fitting model: ",idx
#        print "best misfit                : ",mf
    warg = "-"+ptype+wtype
    p2 = Popen(["gpdcreport",warg,str(mode),"-i",mdlfile],stdin=PIPE,stdout=PIPE).communicate(input=str(idx))[0]
    p = []
    v = []
    for line in p2.split('\n'):
        line=line.rstrip()
        if line.find('#') != -1: continue
        if len(line.split()) < 1: break
        p.append(1./float(line.split()[0]))
        v.append(1./(1000.*float(line.split()[1])))
    return p,v


if __name__ == '__main__':
    usage = "usage: %s filename [options]"%os.path.basename(sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option("-v",dest="verbose",action="store_true",
                      help="print out information on model file to stdout",
                      default=False)
    opts,args = parser.parse_args()
    if len(args) != 1:
        print usage
        sys.exit(1)
    mdlfile = args[0]
    p,v = gpdccurve(mdlfile,wtype="R",ptype="p")
    plot(p,v)
    show()



