#!/usr/bin/env mypython
"""
Create parameter files for the 3D inversion.
"""

import os

dirn = '/data/wanakaII/yannik/cnipse/inversion/3D/paramfiles'
master_tvz = os.path.join(dirn,'Benson_-38.5_176.0.param')
master_ecb = os.path.join(dirn,'Benson_-39.25_177.25.param')
#master_ecb = os.path.join(dirn,'Benson_-39.0_177.0.param')
#master_gw = os.path.join(dirn,'Benson_-38.25_177.0.param')
master_gw = os.path.join(dirn,'Benson_-39.0_176.5.param')

#master_tvz = os.path.join(dirn,'nlayer_model_slow_mantle.param')
#master_ecb = os.path.join(dirn,'nlayer_model_slow_mantle.param')
#master_gw = os.path.join(dirn,'nlayer_model_slow_mantle.param')



### TVZ
tvz = [(175.75,-39.0),(176.0,-39.0),
(175.75,-38.75),(176.0,-38.75),(176.25,-38.75),(176.5,-38.75),
(175.5,-38.5),(175.75,-38.5),(176.0,-38.5),(176.25,-38.5),(176.5,-38.5),
(175.5,-38.25),(175.75,-38.25),(176.0,-38.25),(176.25,-38.25),(176.5,-38.25),(176.75,-38.25),
(175.75,-38.0),(176.0,-38.0),(176.25,-38.0),(176.5,-38.0),(175.5,-39.5)]

eastcoastbasin = [
    (176.5,-39.5),(176.75,-39.5),(177.0,-39.5),(177.25,-39.5),(177.5,-39.5),
    (176.75,-39.25),(177.0,-39.25),(177.25,-39.25),(177.5,-39.25),
    (177.0,-39.0),(177.25,-39.0),(177.5,-39.0),
    (177.0,-38.75),(177.25,-38.75),(177.5,-38.75),
    (177.5,-38.5)]


greywacke = [(175.0, -39.5),(175.0, -39.25),(175.0, -39.0),(175.0, -38.75),(175.0, -38.5),
 (175.0, -38.25),(175.0, -38.0),(175.25, -39.5),(175.25, -39.25),(175.25, -39.0),(175.25, -38.75),
 (175.25, -38.5),(175.25, -38.25),(175.25, -38.0),(175.5, -39.25),(175.5, -39.0),
 (175.5, -38.75),(175.5, -38.0),(175.75, -39.5),(175.75, -39.25),(176.0, -39.5),(176.0, -39.25),
 (176.25, -39.5),(176.25, -39.25),(176.25, -39.0),(176.5, -39.25),(176.5, -39.0),(176.75, -39.0),
 (176.75, -38.75),(176.75, -38.5),(176.75, -38.0),(177.0, -38.5),(177.0, -38.25),(177.0, -38.0),
 (177.25, -38.5),(177.25, -38.25),(177.25, -38.0),(177.5, -38.25),(177.5, -38.0)]

paramdir = '/data/wanakaII/yannik/cnipse/inversion/3D/rayleigh/5-25s/200_100_1/nlayer_slow_mantle_with_error/param'
#paramdir = '/data/wanakaII/yannik/cnipse/inversion/3D/rayleigh/5-25s/200_100_1/3_regions_no_gw_lvz_with_error_30G_1/param'
paramdir = '/data/wanakaII/yannik/cnipse/inversion/3D/love/6-11s/200_100_1/3_regions_no_gw_lvz_with_error_no_u_30G/param'
paramdir = '/data/wanakaII/yannik/cnipse/inversion/3D/joint/200_100_1/3_regions_no_gw_lvz_with_error_no_lc/param'

if not os.path.isdir(paramdir):
    os.makedirs(paramdir)
paramname = '3_regions'
for i in tvz:
    lon,lat = i
    link_nm = os.path.join(paramdir,"%s_%f_%f.param"%(paramname,lat,lon))
    os.symlink(master_tvz,link_nm)

for i in eastcoastbasin:
    lon,lat = i
    link_nm = os.path.join(paramdir,"%s_%f_%f.param"%(paramname,lat,lon))
    os.symlink(master_ecb,link_nm)

for i in greywacke:
    lon,lat = i
    link_nm = os.path.join(paramdir,"%s_%f_%f.param"%(paramname,lat,lon))
    os.symlink(master_gw,link_nm)


#runlat = [-39.5,-39.25,-39.0,-38.75,-38.5,-38.25,-38.0]
#runlon = [175.0,175.25,175.5,175.75,176.0,176.25,176.5,176.75,177.0,177.25,177.5]
#greywacke = []
#for _lon in runlon:
#    for _lat in runlat:
#        i = (_lon,_lat)
#        if i not in tvz and i not in eastcoastbasin:
#            greywacke.append(i)
#for i in tvz:
#    lon,lat = i
#    plot(lon,lat,'ro')
#
#for i in eastcoastbasin:
#    lon,lat = i
#    plot(lon,lat,'bo')
#
#for i in greywacke:
#    lon,lat = i
#    plot(lon,lat,'ko')


