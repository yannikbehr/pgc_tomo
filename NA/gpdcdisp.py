#!/usr/local/bin/python
"""
plot set of dispersion curves color coded by their misfit
"""

import os, os.path, sys, tempfile
from gmtpy import *
from pylab import *

def mkdisp_plot(mdl,disp,ptype='p',wtype='R',mode=0):
    fn = tempfile.mktemp()
    os.system('gpdcreport -%s%s %d %s >%s'%(ptype,wtype,mode,mdl,fn))
    #fn = '/home/behrya/dev/data/tmp.txt'
    nmdl = 1000
    f = open(fn)
    f.readline()
    nemax = int(f.readline().split()[4])
    curpos = f.tell()
    f.readline()
    nrow = 0
    while True:
        line = f.readline()
        if line.find('#') != -1: break
        a = line.split()
        nrow += 1
    f.seek(curpos)
    models = zeros((nemax,nrow,2))
    data   = zeros(nemax)
    ne = 0
    while True:
        line = f.readline()
        if not line: break
        if line.find('value') != -1:
            a = line.split()
            mf = float(a[9].split('=')[1])
            for _i in xrange(nrow):
                line = f.readline()
                a = line.split()
                models[ne,_i,0]=1./float(a[0])
                models[ne,_i,1]=1./(float(a[1])*1000.)
                data[ne]=mf
            ne += 1
    f.close()


    mfmax = data.max()
    mfmin = data.min()
    step = (mfmax-mfmin)/10.
    gmt = GMT()
    rng = '5/30/1.5/5'
    sca = 'X6c/4c'
    anot= 'a5f1/a1f.5WSne'
    scld = '6.5c/0.5c/4c/.25c'
    sclb = '%f:misfit:/::'%(step/.2)
    cptf = tempfile.mktemp('tmp.cpt')
    gmt.makecpt(C='wysiwyg',T='%f/%f/%f'%(mfmin,mfmax,step),I=True,out_filename=cptf)
    fout = 'test.ps'

    fn = tempfile.mktemp()
    #fn = 'tmp.txt'
    #f = open(fn,'w')
    for _i in reversed(data[random_integers(0,nemax-1,nmdl)].argsort()):
    #for _i in reversed(data.argsort()):
    #    f.write('> -Z%f\n'%data[_i])
    #    savetxt(f,models[_i])
        dd = ['> ','-Z %f'%data[_i]]
        dd = reshape(dd,(1,2))
        gmt.psxy(R=rng,J=sca,C=cptf,B=anot,W='2',m=True,in_rows=append(dd,models[_i],axis=0))
    gmt.save(fout)
    os.system('gv %s&'%fout)
    #f.close()
    return fn,data.max(),data.min()

if __name__ == '__main__':
    mdl   = '/home/behrya/dev/data/mt_ray_u_c_reports/mt_ray_u_c_01.report'
    disp  = '/home/behrya/dev/data/dispcurves/MATA_TIKO_ray_c_res.txt'
    fd,mfmax,mfmin = mkdisp_plot(mdl,disp)
#    step = (mfmax-mfmin)/10.
#    os.system('psxy -R%s -J%s -B%s -X.6c -Y.5c -C%s -W2 -m %s -K >%s'%(rng,sca,anot,cptf,fd,fout))
#    os.system('psxy -R -J -B -W3  %s -O -K >>%s'%(disp,fout))
#    os.system('psscale -C%s -D%s -B%s -O >>%s'%(cptf,scld,sclb,fout))
