#!/usr/bin/env python
"""
Plot results from pseudo 3D inversion
"""
from enthought.mayavi.modules.api import Surface, Axes, ScalarCutPlane, Outline, IsoSurface
from enthought.mayavi.sources.api import VTKDataSource
from enthought.tvtk.api import tvtk
from enthought.mayavi import mlab
from scipy.io import savemat,loadmat
from collections import defaultdict
from gmtpy import GMT
import pdb
import sys
import os
import numpy as np
import tempfile
#import pylab as plt

class Plot3DError(Exception): pass

def get_average_mdl(repfile,misfit,sreso,dreso):
    """
    calculate weighted average model for each run
    """
    tmp = tempfile.mktemp()
    os.system('/usr/local/Geopsy.org/bin/gpdcreport '+repfile+' >%s'%tmp)
    models, data, nd, ne, nrow, nlayer = get_models(tmp)
    if misfit <= data.min():
        misfit = median(data)
    print misfit
    print repfile
    ### calculate density
    grid1, grid2, x, y, smean, dmean = dplot(models,data,nd,ne,dreso=dreso, sreso=sreso,mf=misfit,nlayer=nlayer)
    lat = os.path.basename(repfile).split('_')[1]
    lon = os.path.basename(repfile).split('_')[2].split('.report')[0]
    return np.vstack((smean,-dmean)).T

def get_best_model(repfile,dreso,dmax=40.,dmin=0.):
    tmp = tempfile.mktemp()
    ### read in models from dinver-output
    os.system('gpdcreport -best 1 %s >%s'%(repfile,tmp))
    try:
        mdl = np.loadtxt(tmp,skiprows=3)
    except IOError:
        raise Plot3DError("Cannot read report file %s." % repfile)
        return
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
    a = np.array(nmdl)
    a[-1,0] = dmax
    depth = np.arange(dmin,dmax+dreso,dreso)
    nvs = np.interp(depth,a[:,0],a[:,1])
    return np.vstack((nvs,-depth)).T
                
    
def get_1Dresults(runlat,runlon,savedir,wavet,new=False,matfile=None,av=False,
                  best=True,err=False,errdepth=False,group=True,dmax=40.,dreso=0.5,
                  prefix='cnipse'):
    """
    read in average models for each grid-cell
    """
    if new:
        lat = np.array([])
        lon = np.array([])
        depth = np.array([])
        vs = np.array([])
        vp = np.array([])
        vserr = np.array([])
        for _i in xrange(len(runlat)):
            for _j in xrange(len(runlon)):
                if group:
                    fn = '%s/reports/%s_%f_%f.report'%(savedir,prefix,runlat[_i],runlon[_j])
                else:
                    fn = '%s/reports/%s_%f_%f.report'%(savedir,prefix,runlat[_i],runlon[_j])
                print _i,_j,fn
                sreso = 0.05 #horizontal resolution
                dreso = dreso  #vertical resolution
                misfit = 0.2
                if av:
                    prf = get_average_mdl(fn,misfit,sreso,dreso)
                elif best:
                    try:
                        prf = get_best_model(fn,dreso,dmax=dmax,dmin=0.)
                    except Plot3DError:
                        continue
                elif err:
                    prf = get_err_vs_mdl(fn,0.5)
                elif errdepth:
                    prf = get_err_depth(fn,0.5)

                for _d in xrange(prf.shape[0]):
                    lat = np.append(lat,runlat[_i])
                    lon = np.append(lon,runlon[_j])
                    depth = np.append(depth,prf[_d,1])
                    vs = np.append(vs,prf[_d,0])
                    if prf.shape[1] > 2:
                        vserr = np.append(vserr,prf[_d,2])
                        vp = np.append(vp,prf[_d,3])
        savemat(matfile,{'lat':lat,'lon':lon,'depth':depth,'vs':vs,'vserr':vserr,'vp':vp})
        return np.atleast_2d(lat).T,np.atleast_2d(lon).T,np.atleast_2d(depth).T,np.atleast_2d(vs).T,np.atleast_2d(vserr).T,np.atleast_2d(vp).T
    if not new:
        a =  loadmat(matfile,struct_as_record=True)
        lat = a['lat']
        lon = a['lon']
        depth= a['depth']
        vs = a['vs']
        vserr = a['vserr']
        vp = a['vp']
        return lat,lon,depth,vs,vserr,vp

def plot_canada_map():
    """
    Plot coastline as tubes in 3D.
    """
    d = defaultdict(list)
    gmt = GMT(config={'PAGE_ORIENTATION':'landscape'})
    #rng = '-70/-52/40/52.'
    rng = '-150/-49/40./89'
    scl = 'M10c'
    gmt.pscoast(R=rng,J=scl,B='a0.5wsne',D='l',W='thinnest',m=True,A='2/1/1')
    a = gmt.output.getvalue().split('\n')
    z = 0.
    cnt = 0
    connections = list()
    for _l in a:
        if _l.startswith('#'):continue
        if _l.startswith('>'):
            cnt += 1
            continue
        try:
            d[cnt].append(map(float,_l.split('\t')))
        except ValueError:
            print _l

    for _k in d.keys():
        ar = np.array(d[_k])
        x = (6371)*np.cos(2*np.pi*ar[:,1]/360.)*np.cos(2*np.pi*ar[:,0]/360.)
        y = (6371)*np.cos(2*np.pi*ar[:,1]/360.)*np.sin(2*np.pi*ar[:,0]/360.)
        z = (6371)*np.sin(2*np.pi*ar[:,1]/360.)
        pts = mlab.plot3d(x,y,z,tube_radius=2.0,color=(0,0,0))


def get_tvtk_grid(la,lo,dp,vs):
    """
    Convert 1D inversion results into a TVTK structured grid.
    """
    ndepth = np.unique(dp).size
    nla = np.unique(la).size
    nlo = np.unique(lo).size
    points = np.empty((nla*nlo*ndepth,3))
    vals = np.empty(nla*nlo*ndepth)
    # reformat dp,la,lo and vs arrays so that la changes the fastest
    # and dp the slowest (required by mayavi)
    newla = la.reshape(nla,nlo,ndepth)
    newlo = lo.reshape(nla,nlo,ndepth)
    newdp = dp.reshape(nla,nlo,ndepth)
    newvs = vs.reshape(nla,nlo,ndepth)
    cnt = 0
    for _j in xrange(ndepth):
        for _i in xrange(nlo):
            start = cnt*nla
            end = start+nla
            points[start:end,0] = (newdp[:,_i,_j]+6371)*np.cos(2*np.pi*newla[:,_i,_j]/360.)*np.cos(2*np.pi*newlo[:,_i,_j]/360.)
            points[start:end,1] = (newdp[:,_i,_j]+6371)*np.cos(2*np.pi*newla[:,_i,_j]/360.)*np.sin(2*np.pi*newlo[:,_i,_j]/360.)
            points[start:end,2] = (newdp[:,_i,_j]+6371)*np.sin(2*np.pi*newla[:,_i,_j]/360.)
            vals[start:end] = newvs[:,_i,_j]
            cnt += 1
    sgrid = tvtk.StructuredGrid(dimensions=(nla,nlo,ndepth),points=points)
    sgrid.point_data.scalars = vals.copy()
    sgrid.point_data.scalars.name = 'Shear-velocity'
    return sgrid

def convert_pt(lat,lon,depth):
    x = (6371+depth)*np.cos(2*np.pi*lat/360.)*np.cos(2*np.pi*lon/360.)
    y = (6371+depth)*np.cos(2*np.pi*lat/360.)*np.sin(2*np.pi*lon/360.)
    z = (6371+depth)*np.sin(2*np.pi*lat/360.)
    return x,y,z
    

def plot_3d(lats,lons,depth,vs,runlat,runlon,vtkname='rayleigh_c_u.vtk',
            annotate_depth=False,coastline=False,annotate_lat=True,annotate_lon=True):
    """
    plot 3d results using mayavi
    """
    if coastline:
        plot_canada_map()
    sgrid = get_tvtk_grid(lats,lons,depth,vs)
    d = mlab.pipeline.add_dataset(sgrid)
    #sf = mlab.pipeline.surface(d)
    #gx = mlab.pipeline.grid_plane(d)
    #gx.grid_plane.axis = 'x'
    gy = mlab.pipeline.grid_plane(d)
    gy.grid_plane.axis = 'x'
    gy.module_manager.scalar_lut_manager.show_scalar_bar = True
    gy.module_manager.scalar_lut_manager.lut_mode = 'jet'
    gy.module_manager.scalar_lut_manager.data_range = np.array([ 2. ,  4.8])
    gy.module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array([100000, 100000])
    gy.module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array([1, 1])
    gy.module_manager.scalar_lut_manager.scalar_bar_representation.position2 = np.array([ 0.08796009,  0.56264591])
    gy.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.03396896,  0.39182879])
    gy.actor.mapper.progress = 1.0
    gy.actor.mapper.scalar_range = np.array([ 0.,  1.])
    gy.actor.mapper.scalar_visibility = True
    gy.actor.property.representation = 'surface'
    gy.grid_plane.position = 6

    #gz = mlab.pipeline.grid_plane(d)
    #gz.grid_plane.axis = 'z'
    if annotate_lat:
        for lat in runlat:
            x,y,z = convert_pt(lat,-58.,10.)
            txt = mlab.text3d(x,y,z,'%d'%(lat),color=(0,0,0),line_width=10.0)
            txt.scale = [20,20,20]
    if annotate_lon:
        for lon in runlon[1::]:
            x,y,z = convert_pt(49.,lon,10.)
            txt = mlab.text3d(x,y,z,'%d'%(lon),color=(0,0,0),line_width=10.0)
            txt.scale = [20,20,20]
    if annotate_depth:
        for dp in [-10,-40,-80,-120]:
            x,y,z = convert_pt(49,-68.,dp)
            txt = mlab.text3d(x,y,z,'%d km'%(dp),color=(0,0,0),line_width=10.0)
            txt.scale = [20,20,20]

    ### Include 3D screenshot in matplotlib
    #arr = mlab.screenshot()
    #import pylab as pl
    #pl.imshow(arr)
    #pl.show()
    mlab.text(0.76,0.86,'49N',width=0.1)
    mlab.show()

def ascii_output(lon,lat,depth,vs,fout):
    f = open(fout,'w')
    for i in xrange(lon.size):
        print >>f,'%4.3f  %4.3f  %4.3f  %4.3f'%(lon[i],lat[i],depth[i],vs[i])
    f.close()

    
def main(new=True):
    ### St. Lawrence Basin:
    #runlat=np.arange(43.,52.)  
    #runlon=np.arange(-68.,-55.)
    ### Canada:
#    runlat = np.arange(40.,89.)
#    runlon = np.arange(-150.,-49.)
#    savedir = '../../results/canada'
#    result = 'canada_points.npz'
    
    ### Canada 2nd run
    runlat = np.arange(40.,71.)
    runlon = np.arange(-150.,-54.)
    savedir = '../../results/canada_2nd_run'
    result = 'canada_points_2nd_run.npz'
    
    matfile = os.path.join(savedir,'canada.mat')
    if new:
        lat,lon,depth,vs,vserr,vp = get_1Dresults(runlat,runlon,savedir,'rayleigh',new=True,matfile=matfile,
                                                  dmax=120.,dreso=1.0,prefix='canada')
        np.savez(result,lat=lat,lon=lon,depth=depth,vs=vs,vserr=vserr,vp=vp)
    else:
        npz = np.load(result)
        lat = npz['lat']
        lon = npz['lon']
        depth = npz['depth']
        vs = npz['vs']
    #ascii_output(lon,lat,depth,vs,'canada_surface_tomo_2nd_run.txt')
    plot_3d(lat,lon,depth,vs,runlat,runlon,coastline=True,annotate_depth=True,annotate_lon=True,annotate_lat=False)


if __name__ == '__main__':
    #savedir ='/data/wanakaII/yannik/cnipse/inversion/3D/love/5-11s/200_100_1/3_regions_no_gw_lvz_with_error_no_u'
    #savedir ='/data/wanakaII/yannik/cnipse/inversion/3D/rayleigh/5-25s/200_100_1/3_regions_no_gw_lvz_with_error_30G_1'
    #runlat = [-39.5,-39.25,-39.0,-38.75,-38.5,-38.25,-38.0]
    #runlon = [175.0,175.25,175.5,175.75,176.0,176.25,176.5,176.75,177.0,177.25,177.5]
    #matfile = os.path.join(savedir,'cnipse_3d_rayleigh_err.mat')
    #lat,lon,depth,vs,vserr,vp = get_1Dresults(runlat,runlon,savedir,'rayleigh',new=False,matfile=matfile)
    #x,y = latlon2nzmg(lat,lon)
    #plot_nzmap()
    #plot_3d(x,y,depth,vs,os.path.join(savedir,'cni_rayleigh.png'),new=True,
    #        vtkname=os.path.join(savedir,'cni_rayleigh.vtk'))
    #plot_canada_map()
    main(new=False)
