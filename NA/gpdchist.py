#!/usr/local/bin/python
"""
draw histograms for inversion parameters
"""

from pylab import *
import os, os.path, sys
from gpdcreport2mdl import *
from matplotlib import rcParams
import matplotlib.mlab as mlab
rcParams['figure.figsize']=11.69,8.26
rcParams['figure.dpi']=72
rcParams['figure.subplot.wspace']=.4
rcParams['figure.subplot.hspace']=.4
from scipy import optimize

def convmdl(mdlname):
    os.system('gpdcreport '+mdlname+' >tmp.mdl')
    models, data, nd, ne, nl = get_models('tmp.mdl')
    os.remove('tmp.mdl')
    nmdl = zeros((nd,ne+1))
    for _j in range(ne+1):
        for _k in range(nd):
                nmdl[_k,_j]=models[_j*nd+_k]
    return nmdl,nl,nd

def gauss(x,y):
    fitfunc = lambda p,x: (1/sqrt(2*pi*p[0]**2))*exp(-(x-p[1])**2/(2*p[0]**2))
    errfunc = lambda p,x,y: fitfunc(p,x)-y
    gaussian = lambda m,s,x: (1/sqrt(2*pi*s**2))*exp(-(x-m)**2/(2*s**2))
    p0 = [10.,50.]
    p1, success = optimize.leastsq(errfunc,p0[:],args=(x,y))
    sigm,mean = p1
    return sigm,mean,gaussian(mean,sigm,x)

def main(nmdl,nl,nd,mdlname):
    row = 0
    np = nd/nl
    lbl="th0 th1 th2 th3 th4 vp0 vp1 vp2 vp3 vp4 vs0 vs1 vs2 vs3 vs4\
    p0 p1 p2 p3 p4".split()
    for row in range(np):
        for _i in range(nl):
            subplot(np+1,5,1+_i+row*nl)
            n, bins, patches = hist(nmdl[row*nl+_i,:], 100, facecolor='black')
            setp(getp(gca(),"xticklabels"),color="k",fontsize="x-small")
            setp(getp(gca(),"yticklabels"),color="k",fontsize="x-small")
            sigm,mean,gaussian = gauss(bins[:-1],n)
            print sigm, mean
            xlabel(lbl[row*nl+_i])
    ### special
    subplot(np+1,5,_i+row*nl+2)
    n,bins,patches = hist(nmdl[0,:]+nmdl[1,:]+nmdl[2,:]+nmdl[3,:],100,facecolor='black')
    setp(getp(gca(),"xticklabels"),color="k",fontsize="x-small")
    setp(getp(gca(),"yticklabels"),color="k",fontsize="x-small")
    sigm,mean,gaussian = gauss(bins[:-1],n)
    print sigm, mean
    xlabel('Moho depth [km]')
    n,bins,patches = hist(nmdl[0,:]+nmdl[1,:],100,facecolor='black')
    setp(getp(gca(),"xticklabels"),color="k",fontsize="x-small")
    setp(getp(gca(),"yticklabels"),color="k",fontsize="x-small")
    sigm,mean,gaussian = gauss(bins[:-1],n)
    print sigm, mean
    n,bins,patches = hist(nmdl[0,:]+nmdl[1,:]+nmdl[2,:],100,facecolor='black')
    setp(getp(gca(),"xticklabels"),color="k",fontsize="x-small")
    setp(getp(gca(),"yticklabels"),color="k",fontsize="x-small")
    sigm,mean,gaussian = gauss(bins[:-1],n)
    print sigm,mean
    suptitle(mdlname)
    figname = mdlname.split('.report')[0]+'_hist.ps'
    savefig(figname,orientation='landscape',papertype='A4')
    os.system('gv '+figname+'&')


if __name__ == '__main__':
    try:
        mdl = sys.argv[1]
    except:
        print "usage: %s report-file"%os.path.basename(sys.argv[0])
        sys.exit(1)
    nmdl,nl,nd = convmdl(mdl)
    main(nmdl,nl,nd,mdl)
