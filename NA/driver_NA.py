#!/usr/bin/python

"""
@package driver_NA.py
\brief wrapper for NA_sampler and NA_bayes code

Wrapper script to run NA_sampler on evison (GNS) and NA_bayes
locally. ascii model-files created on evison have to be converted into
nad-files using the 'writenad'-python routines. Finally
density plots and 1D marginal PDF are plotted.
"""

import os, os.path, sys, time, shutil

myvars = {'runum':1111,
          'velmodel':'models/northland13',
          'dispcurves':"['obs_disp/tiko_mata_na_r.txt']",
          'naparam':'na.in'}

s = open('data/na_test.py','w')
############## writing summary file and check
############## that required files exist
s.write("""
import os, os.path, time, shutil
fn       = 'summary.'+str(%(runum)d)
f = open(fn,'w')
print >>f, 'runnumber       : str(%(runum)d)'
print >>f, 'date            : '+time.asctime()
for dcurve in %(dispcurves)s:
    print >>f, 'dispersion curve: '+dcurve
print >>f,''
print >>f, 'velocity model: %(velmodel)s'
print >>f, '=================================='
for dcurve in %(dispcurves)s:
    if not os.path.isfile(dcurve):
        print "cannot find dispersion curve file "+dcurve
        return
if not os.path.isfile('%(velmodel)s'):
    print "cannot find velocity model file %(velmodel)s"
    return
if not os.path.isfile('%(naparam)s'):
    print "cannot find na.in file %(naparam)s"
    return
tmp = open('%(velmodel)s','r').readlines()
for line in tmp:
    print >>f, line
print >>f, ''
print >>f, 'na.in file:'
print >>f, '=================================='
tmp = open('%(naparam)s','r').readlines()
for line in tmp:
    print >>f, line
f.close()


############### creating surface.in file
f = open('surface.in','w')
print >>f,'#'
print >>f,'#'
print >>f,'#'
print >>f,os.path.basename('%(velmodel)s')
print >>f,os.path.dirname('%(velmodel)s')
print >>f,len(%(dispcurves)s)
for dcurve in %(dispcurves)s:
    print >>f, dcurve
f.close()

############## creating sobs.d file
shutil.copy('sobs_master.d','sobs.d')
f = open('sobs.d','a')
for dcurve in %(dispcurves)s:
    print >>f, dcurve
f.close()
"""%myvars)
s.close()
