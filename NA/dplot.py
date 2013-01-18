#!/usr/local/bin/python
"""
plot density of output model-ensemble from Neighbourhood
Algorithm on a velocity-depth grid
"""

from ctypes import *
from numpy import *
import sys, os, os.path
import pdb

path = os.path.dirname(__file__)
def dplot(models,data,nd,ne,lib=os.path.join(path,'densityp.so'),dmin=0.0,dmax=40.0,
          dreso=0.1,smin=0.0,smax=5.0,sreso=0.05,mf=0.1,prmidx=2,nlayer=5):
    dp = cdll.LoadLibrary(lib)
    grid1 = zeros((int((dmax-dmin)/dreso)+1,int((smax-smin)/sreso)+1))
    grid2 = zeros((int((dmax-dmin)/dreso)+1,int((smax-smin)/sreso)+1))
    smean = zeros(700)
    wtmean = zeros(700)
    idmax = c_int()
    n = dp.dplot(models.ctypes.data_as(POINTER(c_double)),data.ctypes.data_as(POINTER(c_double)),
                 grid1.ctypes.data_as(POINTER(c_double)),grid2.ctypes.data_as(POINTER(c_double)),
                 c_int(ne),c_int(nd),c_int(nlayer),c_int(prmidx),c_double(dmax),c_double(dmin),c_double(dreso),
                 c_double(smax),c_double(smin),c_double(sreso),c_double(mf),smean.ctypes.data_as(POINTER(c_double)),
                 wtmean.ctypes.data_as(POINTER(c_double)),byref(idmax))
    x = arange(smin,smax,sreso)
    y = arange(dmin, dmax, dreso)
    dmean = y[0:idmax.value]
    smean = smean[0:idmax.value]
    return grid1, grid2, x, y, smean, dmean

def dplotpy(models,data,nd,ne,dmin=0.0,dmax=40.0,
            dreso=0.1,smin=0.0,smax=5.0,sreso=0.05,mf=0.1):
    grid1 = zeros((int((dmax-dmin)/dreso)+1,int((smax-smin)/sreso)+1))
    grid2 = zeros((int((dmax-dmin)/dreso)+1,int((smax-smin)/sreso)+1))
    print grid1.shape
    prmidx = 2
    nlayer = 5
    for _i in xrange(ne):
        dh = 0
        for _j in xrange(nlayer):
            th = models[nd*_i+_j]
            v = models[nd*_i+prmidx*nlayer+_j]
            if v> smax: v = smax
            dold = dh
            dh = dh+th
            if dh > dmax: dh = dmax
            id1 = int((dold-dmin)/dreso)
            id2 = int((dh-dmin)/dreso)
            iv1 = int((v-smin)/sreso)
            #print dold, dh, v
            #print id1,id2,iv1
            if id2 > id1:
                for _k in xrange(id1,id2,1):
                    if data[_i] <= mf: grid1[_k,iv1] +=1
                    grid2[_k,iv1] +=1
            if _j > 0:
                iv1 = int((vold-smin)/sreso)
                iv2 = int((v-smin)/sreso)
                if iv2 > iv1:
                    for _k in xrange(iv1,iv2,1):
                        if data[_i] <= mf: grid1[id1,_k] +=1
                        grid2[id1,_k] +=1
                if iv1 > iv2:
                    print vold, v, dh
                    print iv1,iv2,id2
                    for _k in xrange(iv2,iv1,-1):
                        if data[_i] <= mf: grid1[id1,_k] +=1
                        grid2[id1,_k] +=1
            vold = v
    y = arange(dmin, dmax, dreso)
    x = arange(smin,smax,sreso)
    return grid1, grid2, x, y
                    
if __name__ == '__main__':
    try:
        mdlname = sys.argv[1]
    except:
        print "usage: %s dinver-output"%os.path.basename(sys.argv[0])
        sys.exit(1)
        
    from gpdcreport2mdl import *
    import tempfile
    from subprocess import *
    ### read in models from dinver-output
    os.system('gpdcreport '+mdlname+' >tmp.mdl')
    models, data, nd, ne, nrow = get_models('tmp.mdl')
    ### calculate density
    sreso = 0.05 #horizontal resolution
    dreso = 0.2  #vertical resolution
    misfit = 0.1 #misfit threshhold
    grid1, grid2, x, y, smean, dmean = dplot(models,data,nd,ne,dreso=dreso, sreso=sreso,mf=misfit)
    ### write output into temporary files that can be plotted by gmt
    xyz  = tempfile.mktemp()
    xyz2 = tempfile.mktemp()
    grd  = tempfile.mktemp()
    grdcpt = tempfile.mktemp()
    fout   = 'test.ps'
    f = open(xyz,'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if grid1[ii,jj]>0.0:
                print >>f, x[jj],y[ii], grid1[ii,jj]
    f.close()
    f = open(xyz2,'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if grid2[ii,jj] > 0:
                print >>f, x[jj], y[ii], '0.5'
    f.close()

    ### gmt plot
    rng = '-R1/5/0/40'
    scl = '-JX10.2/-12'
    step='-I%f/%f'%(sreso,dreso)
    os.system('xyz2grd %(xyz)s -G%(grd)s %(rng)s %(step)s'%vars())
    os.system('grd2cpt %(grd)s -Cwysiwyg -Z > %(grdcpt)s'%vars())
    os.system('psmask %(xyz2)s %(rng)s %(scl)s %(step)s -Glightgray -K > %(fout)s'%vars())
    os.system('grdimage %(grd)s %(rng)s %(scl)s -Q -C%(grdcpt)s -K -O >> %(fout)s'%vars())
    os.system('psmask -C -K -O >> %(fout)s'%vars())
    os.system("psbasemap %(rng)s %(scl)s -Ba1f.5:'S-velocity [km/s]':/a10f5:'Depth [km]':WnSe -K -O >> %(fout)s"%vars())
    p = Popen('psxy -R -J -B -W3,black -O >> %(fout)s'%vars(),shell=True,stdin=PIPE).stdin
    ### plot weighted average for all models with misfit better than defined threshhold
    for _v,_d in zip(smean,dmean):
        print >>p,_v,_d
    p.close()
    os.system('gv %(fout)s&'%vars())
