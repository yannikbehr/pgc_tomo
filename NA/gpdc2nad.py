#!/usr/local/bin/python

"""
convert output from M. Wathelet's program 'dinver' to
Sambridge format
"""

import sys, os, os.path
sys.path.append(os.path.join(os.environ['PROC_SRC'],'py_ext/na/writenad'))
from mywritenad import WriteNad, WriteNadError
#sys.path.append('/home/behrya/dev/proc-scripts_git/py_ext/na/density')
from pylab import *
from gmtpy import GMT
from gpdcparam2model import ProFile
from gpdcreport2mdl import *

def gpdc2nad(mdlname,prmfn,fout,ns,ns1,nit):
    os.system('gpdcreport '+mdlname+' >tmp.mdl')
    models, data, nd, ne, nrow = get_models('tmp.mdl')
    os.remove('tmp.mdl')
    prf = ProFile(prmfn)
    ranges = zeros((2,nd))
    ranges[0] = prf.thmin+prf.vpmin+prf.vsmin+prf.rhomin
    ranges[1] = prf.thmax+prf.vpmax+prf.vsmax+prf.rhomax
    scales = zeros(nd+1)
    scales[0] = -1
    for _i in xrange(len(ranges[0])):
        scales[_i+1] = ranges[1,_i] - ranges[0,_i]

    if (ns*nit+ns1) != (ne+1):
        print "number of models and na parameters inconsistent"
    WriteNad(fout,mdls=models,data=data,ranges=ranges,scales=scales,
             ns1=ns1,ns=ns,nr=nit,nd=nd,ne=ne)
    

if __name__ == '__main__':
    prmfn = '/home/behrya/dev/data/models/northland13_3.param'
    mdlname = '/home/behrya/dev/data/env_reports/run_06.report'
    fout = 'models_run_06.nad'
    ns = 100
    ns1 = 50
    nit = 100
    gpdc2nad(mdlname,prmfn,fout,ns,ns1,nit)

