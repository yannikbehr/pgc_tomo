#!/usr/local/bin/python

"""
convert output from M. Wathelet's program 'dinver' to
Sambridge format
"""

import sys, os, os.path
sys.path.append('/home/behrya/dev/proc-scripts_git/py_ext/na/writenad')
sys.path.append('/home/behrya/dev/proc-scripts_git/py_ext/na/density')
from subprocess import *
from gmtpy import GMT
import dplot
from gpdccurve import *
from run_evison import *
from gpdcreport2mdl import *
from pylab import *

#from gpdc2nad import gpdc2nad
def main(mdlname,dispdic,title):
    os.system('gpdcreport '+mdlname+' >tmp.mdl')
    models, data, nd, ne, nrow = get_models('tmp.mdl')
    gmt = GMT(config={'BASEMAP_TYPE':'plain','ANOT_FONT_SIZE':8,
                      'LABEL_FONT_SIZE':10,'COLOR_BACKGROUND':'255/255/255',
                      'COLOR_FOREGROUND':'0/0/0','COLOR_NAN':'255/255/255',
                      'PAGE_ORIENTATION':'landscape',
                      'HEADER_FONT_SIZE':15} )
    xyz=gmt.tempfilename('testxyz.txt')
    xyz2=gmt.tempfilename('testxyz2.txt')
    grd=gmt.tempfilename('tmp.grd')
    grdcpt=gmt.tempfilename('tmp.cpt')
    fileout='dens_test.ps'
    rng='1/5/0/40'
    scl='X4.2/-6'
    dreso = 0.2
    sreso = 0.05
    misfit = 0.1
    #grid1, grid2, x, y, smean, dmean = dplot.dplot(models,data,nd,ne,dreso=dreso, sreso=sreso,mf=misfit)
    grid1, grid2, x, y, smean, dmean = dplot.dplot(models,data,nd,ne,dreso=dreso, sreso=sreso,mf=misfit)
    #grid1, grid2, x, y = dplot.dplotpy(models,data,nd,ne,dreso=dreso, sreso=sreso,mf=misfit)
#    matshow(grid1)
#    show()
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
    anot = int(grid1.max()/1000.)*1000/2.
    tick = anot/2
    gmt.xyz2grd(xyz,G=grd,R=rng,I='%f/%f'%(sreso,dreso),out_discard=True)
    gmt.grd2cpt(grd,C="wysiwyg",Z=True,out_filename=grdcpt)
    gmt.psmask(xyz2,R=rng,T=True,J=scl,I='%f/%f'%(sreso,dreso),G='lightgray')
    gmt.grdimage(grd,J=scl,R=rng,Q=True,C=grdcpt)
    gmt.psbasemap(R=rng,J=scl,B='a1f.5:S-velocity [km/s]:/a10f5:Depth [km]::.%s:WnSe'%title)
    gmt.psxy(R=True,J=True,B=True,W='3,black',in_columns=[smean,dmean])
    f = open('/home/behrya/dev/data/mt_fixed_layers_ray_c_u_mean.txt','w')
    for _p,_v in zip(dmean,smean):
        print >>f,_p,_v
    f.close()
    gmt.psscale(C=grdcpt,D='1.0/1./4c/.4ch',B='a%df%d:No. of models:/::'%(anot,tick))
    ### plot dispersion curves
    gmt.psbasemap(R='5/30/2.0/5.0',J='X4.2/2.5',X='5',B='a1f.5:Period [s]:/a1f.5:Velocity [km/s]:WnSe')
    for _d in dispdic.keys():
        vo = load(dispdic[_d][0])
        p,v = gpdccurve(mdlname,wtype=dispdic[_d][1],ptype=dispdic[_d][2])
        gmt.psxy(R=True,J=True,B=True,W='3,black',in_columns=[p,v])
        gmt.psxy(R=True,J=True,B=True,W='3,red',in_columns=[vo[:,0],vo[:,1]])

    gmt.save(fileout)
    os.system('gv '+fileout+'&')
    


if __name__ == '__main__':
    #try:
    #    mdlname = sys.argv[1]
    #except:
    #    print "usage: %s gp-model-file"%os.path.basename(sys.argv[0])
    #    sys.exit(0)
    flist = sys.stdin.read().split('\n')
    _d = {}
    for line in flist:
        if len(line) < 1: continue
        if line[0] == '#': continue
        if line.find('mdln') != -1:
            mdlname  = line.split('=')[1]
        if line.find('dn1') != -1:
            _d['dn1'] = []
            _d['dn1'].append(line.split('=')[1])
        if line.find('wtype1') != -1:
            _d['dn1'].append(line.split('=')[1])
        if line.find('ptype1') != -1:
            _d['dn1'].append(line.split('=')[1])
        if line.find('dn2') != -1:
            _d['dn2'] = []
            _d['dn2'].append(line.split('=')[1])
        if line.find('wtype2') != -1:
            _d['dn2'].append(line.split('=')[1])
        if line.find('ptype2') != -1:
            _d['dn2'].append(line.split('=')[1])
        if line.find('title') != -1:
            title = line.split('=')[1]
        if line[0]=='>':
            main(mdlname,_d,title)
            _d = {}



