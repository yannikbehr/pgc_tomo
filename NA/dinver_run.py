#!/usr/bin/env mypython
"""
run dinver from script
"""

import os
import sys
from numpy import *
import subprocess as sp
from gpdcreport2mdl import *
import tempfile
sys.path.append(os.path.dirname('__file__'))
from dplot import dplot
from subprocess import *
import glob
import pdb

class DinverException(Exception): pass

def get_best_model(repfile,dreso,dmax=40.,dmin=0.):
    tmp = tempfile.mktemp()
    ### read in models from dinver-output
    os.system('/usr/local/Geopsy.org/bin/gpdcreport -best 1 %s >%s'%(repfile,tmp))
    mdl = loadtxt(tmp,skiprows=3)
    dh = 0
    nmdl = []
    nlayer = mdl.shape[0]
    for i in xrange(nlayer):
        th = mdl[i,0]/1000.
        vs = mdl[i,2]/1000.
        dold = dh
        nmdl.append([dold,vs])
        dh += th
        nmdl.append([dh,vs])
    a = array(nmdl)
    a[-1,0] = dmax
    depth = arange(dmin,dmax+dreso,dreso)
    nvs = interp(depth,a[:,0],a[:,1])
    return vstack((nvs,-depth)).T
                
    
def plot_rep(repfile,pixfile,realmod):
    ### read in models from dinver-output
    os.system('/usr/local/Geopsy.org/bin/gpdcreport '+repfile+' >tmp.mdl')
    models, data, nd, ne, nrow, nlayer = get_models('tmp.mdl')
    #os.remove('tmp.mdl')
    ### calculate density
    sreso = 0.05 #horizontal resolution
    dreso = 0.5  #vertical resolution
    misfit = 1.0 #misfit threshhold
    grid1, grid2, x, y, smean, dmean = dplot(models,data,nd,ne,dreso=dreso,sreso=sreso,mf=misfit,nlayer=nlayer)
    ### write output into temporary files that can be plotted by gmt
    xyz  = tempfile.mktemp()
    xyz  = 'xyz.txt'
    xyz2 = tempfile.mktemp()
    grd  = tempfile.mktemp()
    grdcpt = tempfile.mktemp()
    fout   = pixfile
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
    p = Popen('psxy -R -J -B -W3,black -O -K>> %(fout)s'%vars(),shell=True,stdin=PIPE).stdin
    ### plot weighted average for all models with misfit better than defined threshhold
    for _v,_d in zip(smean,dmean):
        print >>p,_v,_d
    p.close()
    p = Popen('psxy -R -J -B -W3,red -O -K >> %(fout)s'%vars(),shell=True,stdin=PIPE).stdin
    ### plot true model
    for _v,_d in zip(realmod[:,0],realmod[:,1]):
        print >>p,_v,_d
    p.close()
    prf = get_best_model(repfile,dreso)
    p = Popen('psxy -R -J -B -W3,white -O >> %(fout)s'%vars(),shell=True,stdin=PIPE).stdin
    ### plot best model
    for l in prf:
        print >>p,l[0],-l[1]
    p.close()
    os.system('gv %(fout)s&'%vars())



def get_disp(model,fmin,fmax,gv=False,gpdcbin='/usr/local/Geopsy.org/bin/gpdc'):
    """
    calculate synthetic dispersion curve using geopsy's gpdc
    """
    mfile = 'model.txt'
    dfile = 'disp_curve.txt'
    if os.path.isfile(mfile):
        os.remove(mfile)
    if os.path.isfile(dfile):
        os.remove(dfile)
    open(mfile,'w').writelines(model)
    nfreq = 100
    if gv:
        cmd = "%(gpdcbin)s -R 1 -L 1 -min %(fmin)f -max %(fmax)f -group -n %(nfreq)d < %(mfile)s > %(dfile)s"%vars()
    else:
        cmd = "%(gpdcbin)s -R 1 -L 1 -min %(fmin)f -max %(fmax)f -n %(nfreq)d < %(mfile)s > %(dfile)s"%vars()
    os.system(cmd)
    a = loadtxt(dfile)
    rw = a[0:nfreq]
    lw = a[nfreq:]
    os.remove(mfile)
    os.remove(dfile)
    return rw, lw

def make_target(wave,ftype,fout,minmf=0.,gv=False,err=None,new=True):
    """
    create target file needed for dinver
    """
    tmpout = 'inputtarget.txt'
    if os.path.isfile(tmpout):
        os.remove(tmpout)
    if new:
        fl = glob.glob(fout+'*')
        if len(fl)>0:
            for _f in fl:
                os.remove(_f)
    r,c = wave.shape
    if err is None:
        savetxt(tmpout,hstack((wave,zeros((r,1)))))
    else:
        savetxt(tmpout,hstack((wave,err.reshape(r,1))))
    if gv:
        if ftype == 'rayleigh':
            cmd = "/usr/local/Geopsy.org/bin/gptarget A -group -R 0 %(fout)s -min-misfit %(minmf)f< %(tmpout)s"%vars()
            os.system(cmd)
        if ftype == 'love':
            cmd = "/usr/local/Geopsy.org/bin/gptarget A -group -L 0 %(fout)s -min-misfit %(minmf)f< %(tmpout)s"%vars()
            os.system(cmd)
    else:
        if ftype == 'rayleigh':
            cmd = "/usr/local/Geopsy.org/bin/gptarget A -R 0 %(fout)s -min-misfit %(minmf)f< %(tmpout)s"%vars()
            os.system(cmd)
        if ftype == 'love':
            cmd = "/usr/local/Geopsy.org/bin/gptarget A -L 0 %(fout)s -min-misfit %(minmf)f< %(tmpout)s"%vars()
            os.system(cmd)
    os.remove(tmpout)

def run_dinver(target,pfile,repf,nit=50,ns0=50,ns=50,nr=50):
    if os.path.isfile(repf):
        os.remove(repf)
    cmd = "/usr/local/Geopsy.org/bin/dinver -i DispersionCurve -optimization -target %(target)s\
    -param %(pfile)s -itmax %(nit)d -ns0 %(ns0)d -ns %(ns)d -nr %(nr)d -o %(repf)s -f 2>/dev/null"%vars()
    ret = os.system(cmd)
    if ret != 0:
        raise DinverException("dinver run failed")


def add_noise(wave):
    f = wave[:,0]
    s = wave[:,1]
    ns = normal(s,0.00001)
    return vstack((f,ns)).T
    
def mdl2tr(model):
    mdl = model.split('\n')
    nl = int(mdl[1])
    d = [0.0]
    vs = []
    _n = 0
    cnt = 0
    while _n < nl:
        a = mdl[_n+2].split()
        depth = (d[cnt]+float(a[0])/1000.)
        d.append(depth)
        d.append(depth)
        vs.append(float(a[2])/1000.)
        vs.append(float(a[2])/1000.)
        cnt += 2
        _n += 1
    return vstack((vs,d[:-1])).T
    

if __name__ == '__main__':
    ### first create synthetic dispersion curve
    if True:
        ### crust 2.0 model for CNIPSE over a halfspace
        model = """
6
1000.0  2500.0    800.0    2100.0 
1000.0  4000.0    2100.0    2400.0 
16000.0 6000.0    3500.0    2700.0 
8000.0  6600.0    3700.0    2900.0 
9000.0  7200.0    4000.0    3050.0 
0.      8000.0    4700.0    3300.0    
"""
    
    if False:
        ### crust 2.0 model for CNIPSE + upper mantle from PREM
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
    fout = 'crust2.0.target'
    paramfile = 'crust2.0.param'
    result = 'crust2.0.report'
    psfile = 'crust2.0d_ensity_with_noise_more_samples.ps'

    rayc,lovc = get_disp(model,0.05,0.2)
    rayu,lovu = get_disp(model,0.05,0.2,gv=True)
    savetxt('crust2.0_love_c.txt',(1./lovc)[::-1])
    #plot(1./rayc[:,0],1./rayc[:,1])
    rayc_n = add_noise(rayc)
    rayu_n = add_noise(rayu)
    lovc_n = add_noise(lovc)
    lovu_n = add_noise(lovu)
    #plot(1./rayc_n[:,0],1./rayc_n[:,1])
    if os.path.isfile(fout):
        os.remove(fout)
    #make_target(rayc_n,'rayleigh',fout)
    #make_target(rayu_n,'rayleigh',fout,gv=True)
    make_target(lovc,'love',fout,minmf=0.1)
    make_target(lovu,'love',fout,gv=True,new=False,minmf=0.1)
    run_dinver(fout,paramfile,result,nit=100,ns0=50,ns=100,nr=50)
    plot_rep(result,psfile,mdl2tr(model))

