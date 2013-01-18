#!/usr/bin/env python
"""
Measure dispersion curves at each geographical point in the surface
wave map and run dinver to find corresponding 1D S-velocity profile at
each point.
"""
from numpy import *
from collections import defaultdict
from scipy.signal import convolve2d
from scipy.io import savemat,loadmat
import scipy.interpolate as scint
import sys
import os
from subprocess import *
sys.path.append('../NA')
sys.path.append('../')
import dinver_run as dr
import param_canada
import pdb
from subprocess import *
#from matplotlib import rcParams
#rcParams = {'backend':'Agg'}
from plot_1D import plot_rep
#from plot_3D import plot_3d, get_1Dresults, latlon2nzmg

DEBUG = True

def plot_disps(lt,ln,lats,lons,smoothmaps,allmaps,periods):
    """
    Plot averaged dispersion curves and dispersion curves that
    contributed to the average.
    """
    ilat = lats.index(lt)
    ilon = lons.index(ln)
    disp = smoothmaps[:,ilat,ilon]
    _p = logspace(log10(min(periods)),log10(max(periods)))
    rep = scint.splrep(periods,disp)
    #ndisp = interp(_p,periods,disp)
    ndisp = scint.splev(_p,rep)
    plot(_p,ndisp,'bo',label='-39.0, 176.0')
    idx = [(ilon,ilat+1),(ilon-1,ilat),(ilon+1,ilat),(ilon,ilat-1)]
    for _i,_j in idx:
        plot(periods,allmaps[:,_i,_j],'k')
    plot(periods,allmaps[:,ilat,ilon],'r')
    xlabel('Period [s]')
    ylabel('Velocity [km/s]')
    #savefig('comp_disp.pdf')

def plot_maps_3d(lons,lats,allmaps):
    """
    Plot surface wave dispersion curve maps at discrete periods as
    a cube of iso-surfaces
    """
    nlon,nlat = meshgrid(lons,lats)
    p = ones((allmaps.shape))
    la = ones((allmaps.shape))
    lo = ones((allmaps.shape))
    cnt = 4.
    for _i in xrange(p.shape[0]):
        la[_i] = nlat
        lo[_i] = nlon
        p[_i] = p[_i]*cnt
        cnt +=1

    src = mlab.pipeline.scalar_field(p,la,lo,allmaps)
    mlab.pipeline.surface(src)
    mlab.axes(xlabel='Period [s]',ylabel='Latitude ',zlabel='Longitude ',
              nb_labels=4)



def load_maps(filenames,new=False,matfile=None,usecols=(0,1,2)):
    """
    load surface wave velocity maps at discrete periods into cube
    """
    d = defaultdict(lambda: defaultdict(dict))
    e = defaultdict(lambda: defaultdict(dict))
    nd = defaultdict(lambda: defaultdict(dict))
    if new:
        mp = loadtxt(filenames[0],usecols=usecols)
        for l in mp[:]:
            d[l[0]][l[1]] = l[2]
            e[l[1]] = 0.
        smoothmaps = zeros((len(filenames),len(e.keys()),len(d.keys())))
        allmaps = zeros((len(filenames),len(e.keys()),len(d.keys())))

        for _j in xrange(len(filenames)):
            ### read in maps into dictionaries
            mp = loadtxt(filenames[_j])
            v = array([])
            for l in mp[:]:
                d[l[0]][l[1]] = l[2]
                e[l[1]] = 0.
                v = append(v,l[2])

            ### reformat dictionaries into matrices
            a = ones((len(e.keys()),len(d.keys())))*v.mean()
            lats = e.keys()
            lons = d.keys()
            lats.sort()
            lons.sort()
            for _i in xrange(len(lats)):
                lat = lats[_i]
                for _k in xrange(len(lons)):
                    lon = lons[_k]
                    if lat not in d[lon].keys():
                        a[_i,_k] = v.mean()
                    else:
                        a[_i,_k] = d[lon][lat]
            ### convolve with smoothing kernel
            #kern = ones((3,3))*1./9.
            kern = array([[0,1,0],[1,2,1],[0,1,0]])/6.
            b = convolve2d(a,kern,'same',fillvalue=v.mean())
            ### save as .mat files
            smoothmaps[_j] = b
            allmaps[_j] = a
            if matfile is not None:
                savemat(matfile,{'smoothmaps':smoothmaps,'allmaps':allmaps,'lats':lats,'lons':lons})
    if not new:
        if matfile is not None:
            allmaps = loadmat(matfile,struct_as_record=True)['allmaps']
            smoothmaps = loadmat(matfile,struct_as_record=True)['smoothmaps']
            lats = loadmat(matfile,struct_as_record=True)['lats'].tolist()
            lons = loadmat(matfile,struct_as_record=True)['lons'].tolist()
        else:
            print "not .mat file given"
            return 0
    return smoothmaps,allmaps,lats,lons


def invert_1d(runlon,runlat,lats,lons,periodsr,periodsl,savedir,paramdir,paramname,
              mapsc=None,mapserrc=None,mapsu=None,mapserru=None,mapslc=None,mapserrlc=None,
              mapslu=None,mapserrlu=None,dryrun=False,plot=True,wavet='rayleigh',misfit=0.1,
              mkparam=None):
    """
    Invert dispersion curves at given latitudes and longitudes for 1D
    S-velocity profiles.
    """
    if not DEBUG:
        widgets = ['dinver: ', pg.Percentage(), ' ', pg.Bar('#'),
                   ' ', pg.ETA()]
        pbar = pg.ProgressBar(widgets=widgets, maxval=len(runlat)*len(runlon)).start()
    if not os.path.isdir(savedir):
        os.makedirs(savedir)
    if not os.path.isdir(os.path.join(savedir,'targets')):
        os.makedirs(os.path.join(savedir,'targets'))
    if not os.path.isdir(os.path.join(savedir,'reports')):
        os.makedirs(os.path.join(savedir,'reports'))
    if not os.path.isdir(os.path.join(savedir,'plots')):
        os.makedirs(os.path.join(savedir,'plots'))
    if not os.path.isdir(paramdir):
        os.makedirs(paramdir)
    
    cnt = 0
    for _ln in runlon:
        for _lt in runlat:
            if not DEBUG:
                pbar.update(cnt)
                cnt += 1
            else:
                print _ln, _lt, "%d/%d"%(cnt+1,len(runlat)*len(runlon))
                cnt += 1
            ilat = lats.index(_lt)
            ilon = lons.index(_ln)
            paramfile = "%s/%s_%f_%f.param"%(paramdir,paramname,lats[ilat],lons[ilon])
            if mkparam is not None:
                mkparam(lats[ilat],lons[ilon],paramfile,lith50=True)
            if not os.path.isfile(paramfile):
                pass
                #print paramfile," does not exist"
                #return
            fout = "%s/targets/canada_%f_%f.target"%(savedir,lats[ilat],lons[ilon])
            result = '%s/reports/canada_%f_%f.report'%(savedir,lats[ilat],lons[ilon])
            psfile = '%s/plots/canada_%f_%f.eps'%(savedir,lats[ilat],lons[ilon])
            if mapsc is not None:
                wavet = 'rayleigh'
                disp = mapsc[:,ilat,ilon]
                _p = logspace(log10(min(periodsr)),log10(max(periodsr)))
                ndisp = interp(_p,periodsr,disp)*1000
                rayc = vstack((1./_p,1./(ndisp))).T
                if mapserrc is not None:
                    errc = mapserrc[:,ilat,ilon]
                    nerrc = interp(_p,periodsr,errc)*1000
                    ### dinver needs the error as slowness
                    ### on the geopsy wiki I found the following formula to conver
                    ### velocity error into slowness error:
                    ### ((1/(V-dv)-1/V)+(1/V-1/(V+dV)))/2
                    nerrc = ((1./(ndisp-nerrc)-1/ndisp)+(1/ndisp-1/(ndisp+nerrc)))/2.
                    dr.make_target(rayc,wavet,fout,minmf=0.01,err=nerrc)
                else:
                    dr.make_target(rayc,wavet,fout,minmf=0.01)
            if mapsu is not None:
                dispu = mapsu[:,ilat,ilon]
                _p = logspace(log10(min(periodsr)),log10(max(periodsr)))
                ndispu = interp(_p,periodsr,dispu)*1000
                rayu = vstack((1./_p,1./(ndispu))).T
                if mapserru is not None:
                    erru = mapserru[:,ilat,ilon]
                    nerru = interp(_p,periodsr,erru)*1000
                    ### dinver needs the error as slowness
                    ### on the geopsy wiki I found the following formula to conver
                    ### velocity error into slowness error:
                    ### ((1/(V-dv)-1/V)+(1/V-1/(V+dV)))/2
                    nerru = ((1./(ndispu-nerru)-1/ndispu)+(1/ndispu-1/(ndispu+nerru)))/2.
                    if mapsc is not None:
                        dr.make_target(rayu,wavet,fout,minmf=0.01,gv=True,new=False,err=nerru)
                    else:
                        dr.make_target(rayu,wavet,fout,minmf=0.01,gv=True,err=nerru)
                else:
                    if mapsc is not None:
                        dr.make_target(rayu,wavet,fout,minmf=0.01,gv=True,new=False)
                    else:
                        dr.make_target(rayu,wavet,fout,minmf=0.01,gv=True)
            if mapslc is not None:
                wavet = 'love'
                disp = mapslc[:,ilat,ilon]
                _p = logspace(log10(min(periodsl)),log10(max(periodsl)))
                ndisp = interp(_p,periodsl,disp)*1000
                lovc = vstack((1./_p,1./(ndisp))).T
                if mapserrlc is not None:
                    errc = mapserrlc[:,ilat,ilon]
                    nerrc = interp(_p,periodsl,errc)*1000
                    ### dinver needs the error as slowness
                    ### on the geopsy wiki I found the following formula to conver
                    ### velocity error into slowness error:
                    ### ((1/(V-dv)-1/V)+(1/V-1/(V+dV)))/2
                    nerrc = ((1./(ndisp-nerrc)-1/ndisp)+(1/ndisp-1/(ndisp+nerrc)))/2.
                    if mapsc is None and mapsu is None:
                        dr.make_target(lovc,wavet,fout,minmf=0.01,err=nerrc)
                    else:
                        dr.make_target(lovc,wavet,fout,minmf=0.01,err=nerrc,new=False)
                else:
                    if mapsc is None and mapsu is None:
                        dr.make_target(lovc,wavet,fout,minmf=0.01)
                    else:
                        dr.make_target(lovc,wavet,fout,minmf=0.01,new=False)
                        
            if mapslu is not None:
                wavet = 'love'
                dispu = mapslu[:,ilat,ilon]
                _p = logspace(log10(min(periodsl)),log10(max(periodsl)))
                ndispu = interp(_p,periodsl,dispu)*1000
                lovu = vstack((1./_p,1./(ndispu))).T
                if mapserrlu is not None:
                    erru = mapserrlu[:,ilat,ilon]
                    nerru = interp(_p,periodsl,erru)*1000
                    ### dinver needs the error as slowness
                    ### on the geopsy wiki I found the following formula to conver
                    ### velocity error into slowness error:
                    ### ((1/(V-dv)-1/V)+(1/V-1/(V+dV)))/2
                    nerru = ((1./(ndispu-nerru)-1/ndispu)+(1/ndispu-1/(ndispu+nerru)))/2.
                    if mapsc is None and mapsu is None and mapslc is None:
                        dr.make_target(lovu,wavet,fout,minmf=0.01,gv=True,err=nerru)
                    else:
                        dr.make_target(lovu,wavet,fout,minmf=0.01,gv=True,err=nerru,new=False)
                else:
                    if mapsc is None and mapsu is None and mapslc is None:
                        dr.make_target(lovu,wavet,fout,minmf=0.01,gv=True)
                    else:
                        dr.make_target(lovu,wavet,fout,minmf=0.01,gv=True,new=False)
                
            if mapsc is None and mapsu is None and mapslc is None and mapslu is None:
                print "No dispersion curves defined"
                return
            if not dryrun:
                if DEBUG:
                    pass
                    #    print paramfile, fout, result
                try:
                    dr.run_dinver(fout,paramfile,result,nit=300,ns0=100,ns=100,nr=100)
                except:
                    print "dinver run failed for ",lats[ilat],lons[ilon]
            if plot:
                try:
                    plot_rep(result,wavet,psfile,show=False,misfit=misfit)
                except Exception, e:
                    print "cannot plot ",result
                    print e
    if not DEBUG:
        pbar.finish()

def get_maps(dirn,alpha,sigma,beta,name,periods):
    """
    get a list of velocity map locations
    """
    f = []
    for _p in periods:
        fn = os.path.join(dirn,'%s_%d_%d_%d.1'%(name,alpha,sigma,_p))
        print fn
        if os.path.isfile(fn):
            f.append(fn)
        else:
            print fn," does not exist"
    return f

def get_res_maps(dirn,alpha,sigma,beta,name,periods):
    """
    get a list of resolution map locations
    """
    f = []
    for _p in periods:
        fn = os.path.join(dirn,'%.1f'%_p,'%d_%d_%d'%(alpha,sigma,beta),'%s_%d.rea'%(name,_p))
        if os.path.isfile(fn):
            f.append(fn)
        else:
            print fn," does not exist"
    return f


if __name__ == '__main__':
    from ConfigParser import SafeConfigParser
    conf = SafeConfigParser()
    try:
        conf.read(sys.argv[1])
    except IndexError:
        print 'usage: %s config_file'%sys.argv[0]
        sys.exit(1)
    runlat = eval(conf.get('area','runlat'))
    runlon = eval(conf.get('area','runlon'))
    savedir = conf.get('data','results')
    dirn = conf.get('data','data')
    alpha = int(conf.get('data','alpha'))
    sigma = int(conf.get('data','sigma'))
    periodsr = eval(conf.get('data','periods'))
    dryrun = conf.getboolean('run-control','testrun')
    doplot = conf.getboolean('run-control','plot')

    #runlat = [ 24.,  25.,  26.,  27.,  28.,  29.,  30.,  31.,  32.,  33.,  34.,
    #           35.,  36.,  37.,  38.,  39.,  40.,  41.,  42.,  43.,  44.,  45.,
    #           46.,  47.,  48.,  49.,  50.,  51.,  52.,  53.,  54.,  55.,  56.,
    #           57.,  58.,  59.,  60.,  61.,  62.,  63.,  64.,  65.,  66.,  67.,
    #           68.,  69.,  70.,  71.,  72.,  73.,  74.,  75.,  76.,  77.,  78.,
    #           79.,  80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,  88.]
    #
    #runlon = [-156., -155., -154., -153., -152., -151., -150., -149., -148.,
    #          -147., -146., -145., -144., -143., -142., -141., -140., -139.,
    #          -138., -137., -136., -135., -134., -133., -132., -131., -130.,
    #          -129., -128., -127., -126., -125., -124., -123., -122., -121.,
    #          -120., -119., -118., -117., -116., -115., -114., -113., -112.,
    #          -111., -110., -109., -108., -107., -106., -105., -104., -103.,
    #          -102., -101., -100.,  -99.,  -98.,  -97.,  -96.,  -95.,  -94.,
    #          -93.,  -92.,  -91.,  -90.,  -89.,  -88.,  -87.,  -86.,  -85.,
    #          -84.,  -83.,  -82.,  -81.,  -80.,  -79.,  -78.,  -77.,  -76.,
    #          -75.,  -74.,  -73.,  -72.,  -71.,  -70.,  -69.,  -68.,  -67.,
    #          -66.,  -65.,  -64.,  -63.,  -62.,  -61.,  -60.,  -59.,  -58.,
    #          -57.,  -56.,  -55.,  -54.,  -53.,  -52.,  -51.,  -50.,  -49.,
    #          -48.,  -47.,  -46.,  -45.,  -44.]
    #runlat = [45.]
    #runlon = [-75.]
    #### St. Lawrence basin (43-51N, 68W-56W)
    #runlat = [43., 44., 45., 46., 47., 48., 49., 50., 51.]
    #runlon = [-68., -67., -66.,  -65.,  -64.,  -63.,  -62.,
    #          -61., -60., -59., -58.,-57., -56.]
    #savedir = '../../results/st_lawrence_basin_canada/'
    #alpha = 1200
    #sigma = 100
    #periodsr = r_[5,6,7,8,9,10,12,15,18,arange(20,95,5)]
    periodsl = r_[5,6,7,8,9,10,12,15,18,arange(20,95,5)]
    if 1:
        if DEBUG:
            print "reading in velocity maps"
        if 1:
            ### Rayleigh waves
            #dirn = '../../data/'
            fnc = get_maps(dirn+'phase/',alpha,sigma,1,'CA',periodsr)
            fnu = get_maps(dirn+'group/',alpha,sigma,1,'CA',periodsr)
            #frc = get_res_maps(dirn+'phase_maps/zz/',alpha,sigma,1,'2lambda_7',periodsr)
            #fru = get_res_maps(dirn+'group_maps/zz/',alpha,sigma,1,'2lambda_7_gv',periodsr)
            #ferrc = get_maps(dirn+'phase_maps_error/zz/',alpha,sigma,1,'2lambda_7_err',periodsr)
            #ferru = get_maps(dirn+'group_maps_error/zz/',alpha,sigma,1,'2lambda_7_gv_err',periodsr)
            wavet = 'rayleigh'
            smoothmapsrc,allmapsrc,lats,lons = load_maps(fnc,new=True)
            smoothmapsru,allmapsru,lats,lons = load_maps(fnu,new=True)
            #smoothmapsrc_err,allmapsrc_err,lats,lons = load_maps(ferrc,new=True)
            #smoothmapsru_err,allmapsru_err,lats,lons = load_maps(ferru,new=True)
            #plot_disps(runlat[0],runlon[0],lats,lons,smoothmapsc,allmapsc,periods)
        if 0:
            ### Love waves
            dirn = '/data/wanakaII/yannik/cnipse/inversion/'
            fnc = get_maps(dirn+'phase_maps/tt/',alpha,sigma,1,'2lambda_5',periodsl)
            fnu = get_maps(dirn+'group_maps/tt/',alpha,sigma,1,'2lambda_5_gv',periodsl)
            frc = get_res_maps(dirn+'phase_maps/tt/',alpha,sigma,1,'2lambda_5',periodsl)
            fru = get_res_maps(dirn+'group_maps/tt/',alpha,sigma,1,'2lambda_5_gv',periodsl)
            ferrc = get_maps(dirn+'phase_maps_error/tt/',alpha,sigma,1,'2lambda_5_err',periodsl)
            ferru = get_maps(dirn+'group_maps_error/tt/',alpha,sigma,1,'2lambda_5_gv_err',periodsl)
            wavet = 'love'
            smoothmapslc,allmapslc,lats,lons = load_maps(fnc,new=True)
            smoothmapslu,allmapslu,lats,lons = load_maps(fnu,new=True)
            smoothmapslc_err,allmapslc_err,lats,lons = load_maps(ferrc,new=True)
            smoothmapslu_err,allmapslu_err,lats,lons = load_maps(ferru,new=True)
            #plot_disps(runlat[0],runlon[0],lats,lons,smoothmapsc,allmapsc,periods)
    if 1:
        if DEBUG:
            print "running inversion"
        paramdir = os.path.join(savedir,'param')
        paramname = 'Lith50_Canada'
        invert_1d(runlon,runlat,lats,lons,periodsr,periodsl,savedir,paramdir,paramname,
                  mapsc=smoothmapsrc,mapsu=smoothmapsru,mapserrc=None,
                  mapserru=None,mapslc=None,mapslu=None,
                  mapserrlc=None,mapserrlu=None,
                  dryrun=dryrun,plot=doplot,wavet=wavet,misfit='all',mkparam=param_canada.crust2mdl)
    if 0:
        if DEBUG:
            print "plotting"
        matfile = os.path.join(savedir,'cnipse_3d_rayleigh.mat')
        lat,lon,depth,vs = get_1Dresults(runlat,runlon,savedir,'rayleigh',new=False,matfile=matfile)
        x,y = latlon2nzmg(lat,lon)
        plot_3d(x,y,depth,vs,os.path.join(savedir,'cni.png'),
                new=False,
                vtkname=os.path.join(savedir,'cni.vtk'))

