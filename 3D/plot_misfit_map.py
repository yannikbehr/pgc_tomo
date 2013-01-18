#!/usr/bin/env mypython
"""
Plot a map of best misfits from 3D inversion runs
"""

from pylab import *
from gmtpy import *
from tempfile import mktemp
import os
import re

DEBUG = True

def get_mf(dirn,runlat,runlon):
    vals = []
    for _i in xrange(len(runlat)):
        for _j in xrange(len(runlon)):
            if DEBUG:
                print runlat[_i],runlon[_j]
            fn = '%s/reports/cnipse_%f_%f.report'%(dirn,runlat[_i],runlon[_j])
            if not os.path.isfile(fn):
                print fn, " does not exist"
                continue
            if DEBUG:
                print fn
            tmp = tempfile.mktemp()
            os.system('gpdcreport -best 1 '+fn+' >%s'%tmp)
            f = open(tmp)
            f.readline()
            l = f.readline()
            match = re.search('\w* value=(-?\d*\.\d*)',l)
            if match:
                mf =  float(match.group(1))
                vals.append([runlon[_j],runlat[_i],mf])
    return vals

if __name__ == '__main__':
    savedir = '../../results/st_lawrence_basin_canada/'
    runlat = [45.0]
    runlon = [-75.0]
    runlat=arange(43.,52.)  
    runlon=arange(-68.,-55.)  
    #vals = get_mf(savedir,runlat,runlon)
    gmt = GMT(config={'HEADER_FONT_SIZE':'24p',
                      'COLOR_NAN':'white'})
    grdf = mktemp()
    grdcpt = mktemp()
    rng = '-70/-50/40/60'
    anot = 'a4f1WnSe:.%s:'%os.path.basename(savedir)
    gmt.xyz2grd(G=grdf,I='1./1.',R=rng,in_rows=vals,out_discard=True)
    gmt.grd2cpt(grdf,E=50,L='0/0.1',C="seis",D=True,out_filename=grdcpt)
    gmt.grdimage(grdf,R=True,J='M8c',C=grdcpt,Y='4c',X='5c')
    #gmt.psbasemap(R=True,J=True,B=anot)
    gmt.psbasemap(R=True,J=True,B=anot)
    gmt.pscoast(R=True,J=True,D='i',W='thinnest')
    gmt.psscale(C=grdcpt,D='3.c/-1.6c/6c/.2ch',B='%f::/::'%.02)
    fout = os.path.join(savedir,"misfit_map.eps")
    foutpdf = fout.replace('.eps','.pdf')
    gmt.save(fout)
    os.system('epstopdf --outfile=%s %s'%(foutpdf,fout))
    os.system("pdfcrop -margins '0. 0. 0. 0.' %s %s "%(foutpdf,foutpdf.replace('.pdf','_crop.pdf')))
    os.system('gv %s'%foutpdf)

