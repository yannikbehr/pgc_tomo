#!/usr/bin/env mypython

import os
import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
from pylab import *
import tempfile
from scipy.io import *
from collections import defaultdict
from gmtpy import GMT
sys.path.append('../NA')


def get_best_model(repfile,dreso,dmax=40.,dmin=0.):
    tmp = tempfile.mktemp()
    ### read in models from dinver-output
    os.system('gpdcreport -best 1 %s >%s'%(repfile,tmp))
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
                
    
def get_1Dresults(runlat,runlon,savedir,wavet,new=False,matfile=None,av=False,
                  best=True,err=False,errdepth=False,group=True):
    """
    read in average models for each grid-cell
    """
    if new:
        lat = array([])
        lon = array([])
        depth = array([])
        vs = array([])
        vp = array([])
        vserr = array([])
        for _i in xrange(len(runlat)):
            for _j in xrange(len(runlon)):
                if group:
                    fn = '%s/reports/cnipse_%f_%f.report'%(savedir,runlat[_i],runlon[_j])
                else:
                    fn = '%s/reports/cnipse_%f_%f.report'%(savedir,runlat[_i],runlon[_j])
                print _i,_j,fn
                sreso = 0.05 #horizontal resolution
                dreso = 0.5  #vertical resolution
                misfit = 0.2
                if av:
                    prf = get_average_mdl(fn,misfit,sreso,dreso)
                elif best:
                    prf = get_best_model(fn,dreso,dmax=40.,dmin=0.)
                elif err:
                    prf = get_err_vs_mdl(fn,0.5)
                elif errdepth:
                    prf = get_err_depth(fn,0.5)

                for _d in xrange(prf.shape[0]):
                    lat = append(lat,runlat[_i])
                    lon = append(lon,runlon[_j])
                    depth = append(depth,prf[_d,1])
                    vs = append(vs,prf[_d,0])
                    if prf.shape[1] > 2:
                        vserr = append(vserr,prf[_d,2])
                        vp = append(vp,prf[_d,3])
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

def latlondep2geo(dp,la,lo,vs):
    ndepth = 81
    nla = 9
    nlo = 13
    #x = zeros((nlo,nla,ndepth))
    #y = zeros((nlo,nla,ndepth))
    #z = zeros((nlo,nla,ndepth))
    #s = zeros((nlo,nla,ndepth))
    points = empty((nla*nlo*ndepth,3))
    vals = empty(nla*nlo*ndepth)
    # reformat dp,la,lo and vs arrays so that la changes the fastest
    # and dp the slowest (required by mayavi)
    #newla = zeros(nla*nlo*ndepth)
    #newlo = zeros(nla*nlo*ndepth)
    #newdp = zeros(nla*nlo*ndepth)
    #newvs = zeros(nla*nlo*ndepth)
    #for _i in xrange(nlo*ndepth):
    #    for _la in xrange(nla):
    #        newla[cnt] = la[_la*nlo*ndepth+_i]
    #        newlo[cnt] = lo[_la*nlo*ndepth+_i]
    #        newdp[cnt] = dp[_la*nlo*ndepth+_i]
    #        newvs[cnt] = vs[_la*nlo*ndepth+_i]
    #        cnt += 1
    newla = la.reshape(nla,nlo,ndepth)
    newlo = lo.reshape(nla,nlo,ndepth)
    newdp = dp.reshape(nla,nlo,ndepth)
    newvs = vs.reshape(nla,nlo,ndepth)
    cnt = 0
    for _j in xrange(ndepth):
        for _i in xrange(nlo):
            start = cnt*nla
            end = start+nla
            points[start:end,0] = (newdp[:,_i,_j]+6371)*cos(2*pi*newla[:,_i,_j]/360.)*cos(2*pi*newlo[:,_i,_j]/360.)
            points[start:end,1] = (newdp[:,_i,_j]+6371)*cos(2*pi*newla[:,_i,_j]/360.)*sin(2*pi*newlo[:,_i,_j]/360.)
            points[start:end,2] = (newdp[:,_i,_j]+6371)*sin(2*pi*newla[:,_i,_j]/360.)
            vals[start:end] = newvs[:,_i,_j]
            cnt += 1
    if 0:
        x = (dp+6371)*cos(2*pi*la/360.)*cos(2*pi*lo/360.)
        y = (dp+6371)*cos(2*pi*la/360.)*sin(2*pi*lo/360.)
        z = (dp+6371)*sin(2*pi*la/360.)
        points[:,1] = y[:,0]
        points[:,2] = x[:,0]
        points[:,0] = z[:,0]
    if 0:
        x = (newdp+6371)*cos(2*pi*newla/360.)*cos(2*pi*newlo/360.)
        y = (newdp+6371)*cos(2*pi*newla/360.)*sin(2*pi*newlo/360.)
        z = (newdp+6371)*sin(2*pi*newla/360.)
        points[:,1] = y[:]
        points[:,0] = x[:]
        points[:,2] = z[:]
    return points,vals


def plot_canada_map():
    d = defaultdict(list)
    gmt = GMT(config={'PAGE_ORIENTATION':'landscape'})
    rng = '-70/-52/40/52.'
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
        ar = array(d[_k])
        x = (6371)*cos(2*pi*ar[:,1]/360.)*cos(2*pi*ar[:,0]/360.)
        y = (6371)*cos(2*pi*ar[:,1]/360.)*sin(2*pi*ar[:,0]/360.)
        z = (6371)*sin(2*pi*ar[:,1]/360.)
        pts = mlab.plot3d(x,y,z,tube_radius=2.0,color=(0,0,0))


def convert_pt(lat,lon,depth):
    x = (6371+depth)*cos(2*pi*lat/360.)*cos(2*pi*lon/360.)
    y = (6371+depth)*cos(2*pi*lat/360.)*sin(2*pi*lon/360.)
    z = (6371+depth)*sin(2*pi*lat/360.)
    return x,y,z
    
if __name__ == '__main__':
    if 1:
        from enthought.tvtk.api import tvtk
        from enthought.mayavi import mlab
        ### St. Lawrence Basin:
        runlat=arange(43.,52.)  
        runlon=arange(-68.,-55.)
        #runlat = [43.,44.,45.]
        #runlon = [-68.,-67.,-66.]
        savedir = '../../results/st_lawrence_basin_canada/'
        matfile = os.path.join(savedir,'st_lawrence_basin.mat')
        #lat,lon,depth,vs,vserr,vp = get_1Dresults(runlat,runlon,savedir,'rayleigh',new=True,matfile=matfile)
        #x,y,z,s = latlondep2geo(depth,lat,lon,vs)
        #np.savez('points.npz',lat=lat,lon=lon,depth=depth,vs=vs,vserr=vserr,vp=vp)
        npz = np.load('points.npz')
        lat = npz['lat']
        lon = npz['lon']
        depth = npz['depth']
        vs = npz['vs']
        pts,s = latlondep2geo(depth,lat,lon,vs)
        sgrid = tvtk.StructuredGrid(dimensions=(9,13,81),points=pts)
        sgrid.points = pts
        sgrid.point_data.scalars = ravel(s.copy())
        sgrid.point_data.scalars.name = 'scalars'
        d = mlab.pipeline.add_dataset(sgrid)
        sf = mlab.pipeline.surface(d)
        gx = mlab.pipeline.grid_plane(d)
        gx.grid_plane.axis = 'x'
        gy = mlab.pipeline.grid_plane(d)
        gy.grid_plane.axis = 'y'
        gz = mlab.pipeline.grid_plane(d)
        gz.grid_plane.axis = 'z'
        for lat in runlat:
            x,y,z = convert_pt(lat,-68.,10.)
            txt = mlab.text3d(x,y,z,'%d'%(lat),color=(0,0,0),line_width=10.0)
            txt.scale = [30,30,30]
        for lon in runlon[1::]:
            x,y,z = convert_pt(43.,lon,10.)
            txt = mlab.text3d(x,y,z,'%d'%(lon),color=(0,0,0),line_width=10.0)
            txt.scale = [30,30,30]
        for dp in [-10,-40]:
            x,y,z = convert_pt(43.,-68.5,dp)
            txt = mlab.text3d(x,y,z,'%d km'%dp,color=(0,0,0),line_width=10.0)
            txt.scale = [30,30,30]
            
        #plot_canada_map()
        mlab.show()
        #arr = mlab.screenshot()
        #import pylab as pl
        #pl.imshow(arr)
        #pl.show()
    
    if 0:
        mpl.rcParams['legend.fontsize'] = 10
        fig = figure()
        ax = Axes3D(fig)
        cset = ax.contour(x[:,:,10], y[:,:,10], z[:,:,10],16)
        cset = ax.contour(x[:,:,20], y[:,:,20], z[:,:,20],16)
    if 0:
        from enthought.tvtk.api import tvtk
        from enthought.mayavi import mlab
        dims = (3,3,81)
        lats = array([43.,44.,])
        lons = array([-68.,-67.,-66.,-65.,-64.])
        depths = linspace(-0,-40,10)
        points = zeros((lats.size*lons.size*depths.size,3))
        #s = empty(lats.size*lons.size*depths.size)
        #s[0:9] = 3.0
        #s[9:18] = 2.0
        s = np.random.random((points.shape[0]))
        cnt = 0
        for d in depths:
            for lat in lats:
                for lon in lons:
                    points[cnt,0] = (d+6371)*cos(pi*lat/180.)*cos(pi*lon/180.)
                    points[cnt,1] = (d+6371)*cos(pi*lat/180.)*sin(pi*lon/180.)
                    points[cnt,2] = (d+6371)*sin(pi*lat/180.)
                    cnt += 1
                    
        sgrid = tvtk.StructuredGrid(dimensions=(lons.size, lats.size, depths.size),points=points)
        #s = np.random.random((dims[0]*dims[1]*dims[2]))
        sgrid.point_data.scalars = ravel(s.copy())
        sgrid.point_data.scalars.name = 'scalars'
        d = mlab.pipeline.add_dataset(sgrid)
        gx = mlab.pipeline.grid_plane(d)
        gy = mlab.pipeline.grid_plane(d)
        gy.grid_plane.axis = 'y'
        gz = mlab.pipeline.grid_plane(d)
        gz.grid_plane.axis = 'z'
        mlab.show()
    if 0:
        from enthought.mayavi import mlab
        mlab.points3d(pts[:,0],pts[:,1],pts[:,2],scale_factor=5.5)
    if 0:
        from enthought.mayavi.scripts import mayavi2
        from enthought.mayavi.sources.vtk_data_source import VTKDataSource
        from enthought.mayavi.modules.api import Outline, GridPlane
        
        mayavi.new_scene()
        src = VTKDataSource(data=sgrid)
        mayavi.add_source(src)
        mayavi.add_module(Outline())
        g = GridPlane()
        g.grid_plane.axis = 'x'
        mayavi.add_module(g)
        g = GridPlane()
        g.grid_plane.axis = 'y'
        mayavi.add_module(g)
        g = GridPlane()
        g.grid_plane.axis = 'z'
        mayavi.add_module(g)

