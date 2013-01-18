#!/usr/bin/env mypython
"""
Plot results from pseudo 3D inversion
"""
from enthought.mayavi import mlab
mlab.options.offscreen = True
from enthought.tvtk.api import tvtk
from enthought.mayavi.sources.api import VTKDataSource
from scipy.io import savemat,loadmat
from collections import defaultdict
from gmtpy import GMT
import pdb
import sys
import os
import numpy as np
sys.path.append('../NA')
from plot_1D import get_best_model
from plot_3D import get_1Dresults, get_tvtk_grid
from gpdcreport2mdl import *
from dplot import dplot
import tempfile
import delaz
import cStringIO
from ConfigParser import SafeConfigParser

       
def convert_pt(lat,lon,depth):
    x = (6371+depth)*np.cos(2*np.pi*lat/360.)*np.cos(2*np.pi*lon/360.)
    y = (6371+depth)*np.cos(2*np.pi*lat/360.)*np.sin(2*np.pi*lon/360.)
    z = (6371+depth)*np.sin(2*np.pi*lat/360.)
    return x,y,z

def get_field(lat,lon,depth,vs,vtkname='vs_model.vtk',visualize=False,new=True):
    ### create velocity field
    if new:
        field = get_tvtk_grid(lat,lon,depth,vs)
        vtksrc = mlab.pipeline.get_vtk_src(field)
        i = VTKDataSource(data=vtksrc[0])
        i.save_output(vtkname)
    else:
        field = mlab.pipeline.open(vtkname,figure=None)
        if visualize:
            edges = mlab.pipeline.extract_edges(field)
            mlab.pipeline.surface(edges, opacity=0.5, line_width=5)
    return field
        
def plot_2d(field,slat,slon,elat,elon,pdepth,new=True):
    """
    Calculate delaunay triangulation and subsequently probe velocity
    field at points of interest.
    """

    ndist = 0.
    nd = 0.
    values = 0.
                    
    ### calculate points between start and end point using a GMT's project program
    fout = 'profile_points.xyp'
    gmt = GMT()
    gmt.project(C='%f/%f'%(slon,slat),E='%f/%f'%(elon,elat),G=50,Q=True,out_filename=fout)
    lon,lat,dist = np.loadtxt(fout,unpack=True)
    os.remove(fout)
    
    cx = []
    cy = []
    cz = []
    cd = []
    cdp = []
    for _lon,_lat,_dist in zip(lon,lat,dist):
        for _d in pdepth:
            x,y,z = convert_pt(_lat,_lon,_d)
            cx.append(x)
            cy.append(y)
            cz.append(z)
            cd.append(_dist)
            cdp.append(_d)
    values = mlab.pipeline.probe_data(field,cx,cy,cz)
    return np.array(cd), np.array(cdp), values

def find_scale(slat,slon,elat,elon):
    gmt = GMT()
    maxdist = 0.
    fout = 'profile_points.xyp'
    for _slat,_slon,_elat,_elon in zip(slat,slon,elat,elon):
        gmt.project(C='%f/%f'%(_slon,_slat),E='%f/%f'%(_elon,_elat),G=100,Q=True,out_filename=fout)
        lon,lat,dist = np.loadtxt(fout,unpack=True)
        if dist.max() > maxdist:
            maxdist = dist.max()
        os.remove(fout)
    return maxdist
    

def test(field):
    savedir = '../../results/canada'
    repfile = os.path.join(savedir,'reports/canada_45.000000_-75.000000.report')
    dreso = 1.0
    bmdl = get_best_model(repfile,dreso,dmax=120.,dmin=0.)
    slat = 50.0
    slon = -75.0
    elat = 50.0
    elon = -120.0
    pdepth = -1*arange(0,120.5,dreso)
    ndist, nd, values = plot_2d(field,slat,slon,elat,elon,pdepth,new=False)
    err = sum((values - bmdl[:,0])**2)/values.size
    
if __name__ == '__main__':
    conf = SafeConfigParser()
    try:
        conf.read(sys.argv[1])
    except IndexError:
        print 'usage: %s config_file'%sys.argv[0]
        sys.exit(1)

    wavet = 'rayleigh'
    runlat = [-39.5,-39.25,-39.0,-38.75,-38.5,-38.25,-38.0]
    runlon = [175.0,175.25,175.5,175.75,176.0,176.25,176.5,176.75,177.0,177.25,177.5]
    savedir = '../../results/canada'
    matfile = os.path.join(savedir,'canada.mat')
    vtkfile = os.path.join(savedir,'canada_3d_%s.vtk'%wavet)
    dreso = 1.0
    pdepth = -1*arange(0,120.5,1.)

    print "reading 3D model" 
    lat,lon,depth,vs,vserr,vp = get_1Dresults(runlat,runlon,savedir,wavet,new=False,matfile=matfile,dmax=120.,dreso=dreso,prefix='canada')
    field = get_field(lat,lon,depth,vs,vtkname=vtkfile,visualize=False,new=True)
#    else:
#        a = np.load(foutprofiles.replace('.eps','.npz'))
#        ndist = a['ndist']
#        nd = a['nd']
#        values = a['values']
#        np.savez(foutprofiles.replace('.eps','.npz'), ndist=ndist, nd=nd, values=values)

    foutprofiles = conf.get('profiles','fout')
    xscale = float(conf.get('profiles','xscale'))
    yscale = float(conf.get('profiles','yscale'))
    direction = conf.get('profiles','direction').lower()
    if direction == 'ns':
        slat = eval(conf.get('ns-profiles','slat'))
        slon = eval(conf.get('ns-profiles','slon'))
        elat = eval(conf.get('ns-profiles','elat'))
        elon = eval(conf.get('ns-profiles','elon'))
        anot = eval(conf.get('ns-profiles','anot'))
    elif direction == 'we':
        slat = eval(conf.get('we-profiles','slat'))
        slon = eval(conf.get('we-profiles','slon'))
        elat = eval(conf.get('we-profiles','elat'))
        elon = eval(conf.get('we-profiles','elon'))
        anot = eval(conf.get('we-profiles','anot'))
    else:
        print "'direction' has to be either 'ns' or 'we'"
        sys.exit(1)
        
    #maxdist = find_scale(slat,slon,elat,elon)
    gmt = GMT(config={'ANOT_FONT_SIZE':14,'LABEL_FONT_SIZE':14,
                      'ANNOT_OFFSET_SECONDARY':'0.1c',
                      'ANNOT_OFFSET_PRIMARY':'0.1c',
                      'LABEL_OFFSET':'0.1c',
                      'FRAME_PEN':'.5p'})
    #scly = (28.-3.-3.*len(slat))/len(slat)
    cptws = gmt.tempfilename('ws.cpt')
    gmt.makecpt(C='seis',D=True,T='3.0/4.5/0.05',out_filename=cptws)
    gmt.psscale(C=cptws,D='8c/.5c/10c/.2ch',B='a%ff.1:Vs:'%(.5))
    cnt = 1
    for _slat,_slon,_elat,_elon in zip(slat,slon,elat,elon):
        print "plotting profile %d"%cnt
        ndist, nd, values = plot_2d(field,_slat,_slon,_elat,_elon,pdepth,new=True)
        lbl1,lbl2 = anot[cnt-1]
        fstr = cStringIO.StringIO()
        for dist,depth,vs in zip(ndist,nd,values):
            fstr.write("%f %f %f\n"%(dist,depth,vs))
        sclx = 18.*ndist.max()/xscale
        scly = yscale
        scl = 'X%fc/%fc'%(sclx,scly)
        print scl
        rng = '%f/%f/%f/%f'%(ndist.min(),ndist.max(),nd.min(),nd.max())
        grdws = gmt.tempfilename('ws.grd')
        grdgd = gmt.tempfilename('grdgd.cpt')
        gmt.xyz2grd(G=grdws,R=rng,I='50./1.',out_discard=True,in_string=fstr.getvalue())
        gmt.grdgradient(grdws,G=grdgd,out_discard=True,A='180',N='e0.8')
        if cnt == 1:
            gmt.grdimage(grdws,I=grdgd,R=True,J=scl,C=cptws,
                         B='a500f100:Distance [km]:/a20f10:Depth [km]:WSne',X='1c',Y='3c')
        else:
            gmt.grdimage(grdws,I=grdgd,R=True,J=scl,C=cptws,
                         B='a500f100:Distance [km]:/a20f10:Depth [km]:WSne',Y='%fc'%(scly+2.5))
        txtstr1 = """%f %f 12 0 1 LB %s"""%(0.0,8.,lbl1)
        txtstr2 = """%f %f 12 0 1 RB %s"""%(ndist.max(),8.,lbl2)
        gmt.pstext(R=True,J=True,G='0/0/0',N=True,in_string=txtstr1)
        gmt.pstext(R=True,J=True,G='0/0/0',N=True,in_string=txtstr2)
        cnt += 1
    gmt.save(foutprofiles)

  
    #test(field)
    
