#!/usr/bin/env mypython
"""
Plot results from pseudo 3D inversion
"""
from scipy.io import savemat,loadmat
from gmtpy import GMT
import sys
import os
import numpy as np
import cStringIO
from ConfigParser import SafeConfigParser

conf = SafeConfigParser()
try:
    conf.read(sys.argv[1])
except IndexError:
    print 'usage: %s config_file'%sys.argv[0]
    sys.exit(1)
ns_slat = eval(conf.get('ns-profiles','slat'))
ns_slon = eval(conf.get('ns-profiles','slon'))
ns_elat = eval(conf.get('ns-profiles','elat'))
ns_elon = eval(conf.get('ns-profiles','elon'))
ns_anot = eval(conf.get('ns-profiles','anot'))

we_slat = eval(conf.get('we-profiles','slat'))
we_slon = eval(conf.get('we-profiles','slon'))
we_elat = eval(conf.get('we-profiles','elat'))
we_elon = eval(conf.get('we-profiles','elon'))
we_anot = eval(conf.get('we-profiles','anot'))
foutmap = conf.get('map','fout')
invres = eval(conf.get('map','invres'))

slat = ns_slat + we_slat
slon = ns_slon + we_slon
elat = ns_elat + we_elat
elon = ns_elon + we_elon
anot = ns_anot + we_anot
# plot cross-section on map
gmt = GMT(config={'ANOT_FONT_SIZE':14,'LABEL_FONT_SIZE':14,
                  'ANNOT_OFFSET_SECONDARY':'0.1c',
                  'ANNOT_OFFSET_PRIMARY':'0.1c',
                  'LABEL_OFFSET':'0.1c',
                  'FRAME_PEN':'.5p',
                  'PLOT_DEGREE_FORMAT':'-D'})
rng = '-145/-50/41/85'
scl = 'L-97.5/63/41/85/15c'
markers='markers.xyp'
gmt.psbasemap(R=rng,J=scl,B='a15f5WSne',Y='2c',X='2c')
gmt.pscoast('-N1',R=True,J=True,D='i',W='0.5p,black',N='2')
cnt = 0
for _slon,_slat,_elon,_elat in zip(slon,slat,elon,elat):
    gmt.project(C='%f/%f'%(_slon,_slat),E='%f/%f'%(_elon,_elat),G=100,Q=True,out_filename=markers)
    gmt.psxy(R=True,J=True,W='2p,red,',in_rows=[[_slon,_slat],[_elon,_elat]])
    #gmt.psxy(markers,R=True,J=True,S='y0.1c',W='1p,blue')
    gmt.psxy(markers,R=True,J=True,S='s0.1c',W='1p,blue')
    lbl1,lbl2 = anot[cnt]
    txtstr1 = """%f %f 10 0 1 RB %s"""%(_slon,_slat,lbl1)
    txtstr2 = """%f %f 10 0 1 LB %s"""%(_elon,_elat,lbl2)
    gmt.pstext(R=True,J=True,W='white',G='blue',N=True,in_string=txtstr1)
    gmt.pstext(R=True,J=True,W='white',G='blue',N=True,in_string=txtstr2)
    cnt += 1
for lon,lat in invres:
    gmt.psxy(R=True,J=True,S='t0.5c',G='red',in_rows=[[lon,lat]])
    
gmt.save(foutmap)
