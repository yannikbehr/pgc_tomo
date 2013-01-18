#!/usr/bin/env mypython
"""
Interpolate the Lith5.0 model onto the grid points of the surface
wave maps.
"""

from gmtpy import GMT
runlat = arange(24.,89.)
runlon = arange(-156.,-43.)
xyfile = 'surface_wave_coord.txt'
xyzfile = 'lith5.0_moho.txt'
f = open(xyfile,'w')
for lat in runlat:
    for lon in runlon:
        print >>f, lon, lat
f.close()

### Plot original
na04_moho = loadtxt('./LITH5.0/NA04_moho.xyf')
gmt = GMT(config={'PAGE_ORIENTATION':'landscape'})
rng = '-145/-50/35/80'
scl = 'L-100/60/45/65/12c'
anot = 'a20f10/a25f5WSne'
fout = 'na04_moho.eps'
tmpgrd = gmt.tempfilename('moho_temp.grd')
mohogrd = gmt.tempfilename('moho.grd')
mohocpt = gmt.tempfilename('moho.cpt')
gmt.xyz2grd(G=mohogrd,I='%f/%f'%(0.25,0.25),R=rng,out_discard=True,in_rows=na04_moho)
gmt.grd2cpt(mohogrd,E=50,L='10/60',C="seis",out_filename=mohocpt)
gmt.grdtrack(xyfile,G=mohogrd,R=True,out_filename=xyzfile)
gmt.grdimage(mohogrd,R=True,J=scl,C=mohocpt)
gmt.pscoast(R=True,J=scl,B=anot,D='i',W='thinnest' )


### Plot resampled 
na04_moho = loadtxt(xyzfile)
gmt.xyz2grd(G=mohogrd,I='%f/%f'%(1.,1.),R=rng,out_discard=True,in_rows=na04_moho)
gmt.grdimage(mohogrd,R=True,J=scl,C=mohocpt,X='13c')
gmt.pscoast(R=True,J=scl,B=anot,D='i',W='thinnest' )
gmt.psscale(C=mohocpt,V=True,D='0c/-.7c/6c/.2ch',B='10::/:km:')
gmt.save(fout)
os.system('gv %s&'%fout)
