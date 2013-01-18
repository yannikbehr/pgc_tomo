#!/usr/bin/env python
"""script to create file from NA results that can be plotted using GMT"""
import numpy, re, string, sys
import pylab as p
import math

class DensPlotGrd:


    def __init__(self):
        self.xgrid = numpy.arange(0,5.01,.01)
        self.ygrid = numpy.arange(0,60,.1)
        self.cols, = self.xgrid.shape
        self.rows, = self.ygrid.shape
        self.grid = numpy.zeros((self.rows,self.cols))
        self.misfitval = []
        self.thickness = []
        self.depth = []
        self.velocity = []
        self.speed = []
        self.sortedmisfit = []
        self.nsublayer = []
        self.ntot = 0
        self.best = 1000
        self.dmin = 0.0
        self.dmax = 40.
        self.dy = .5
        self.vmin = 0.
        self.vmax = 5.
        self.dx = .02
        
    def grd2gmt_alt(self):
        """write out matrix values in columns in order to pipe to gmt"""
        for i in range(0,self.rows):
            for j in range(0,self.cols):
                if self.grid[i][j] == 0:
                    print self.xgrid[j],self.ygrid[i],'NaN'
                else:
                    print self.xgrid[j],self.ygrid[i],self.grid[i][j]
                    

    def val2grd(self, velocity, thickness):
        """look for a value in vel and depth for the corresponding
        bin on the grid spanned by the arrays ygrid and xgrid and count
        up value of these grid points as well as all grid points lying in
        between by 1; result is a 2D matrix where each value corresponds
        to the number of profiles that crosses the corresponding grid-cell"""
        for i in range(0,self.ntot):
            depth = 0
            xold = 0
            yold = 0
            if self.misfitval[i] < self.sortedmisfit[self.best]:
                for k,j  in zip(velocity[i], thickness[i]):
                    depth = depth + j
                    xbin = numpy.searchsorted(self.xgrid,k)
                    ybin = numpy.searchsorted(self.ygrid,depth)
                    #print k,j,xbin,ybin,self.xgrid[xbin],self.ygrid[ybin]
                    if xold != 0 and yold != 0:
                        for l in range(yold, ybin+1):
                            self.grid[l][xold] = self.grid[l][xold]+1
                        if xold > xbin:
                            for m in range(xold-1, xbin, -1):
                                self.grid[l][m] = self.grid[l][m]+1
                        elif xbin > xold:
                            for m in range(xold+1, xbin):
                                self.grid[l][m] = self.grid[l][m]+1
                        else:
                            pass
        
                    xold = xbin
                    yold = ybin



    def get_models(self,filename):
        """get s-velocity structure models from output of
        NA-code"""
        f = open(filename,'r')
        header1 = f.readline()
        header2 = f.readline()
        header3 = f.readline()
        intg = r'\d+'
        real_dn = r'(\d+\.\d*|\d*\.\d+)'
        anywhite = r'\s+'
        anychar = r'\.*'
        pattern1 = r'iteration:\s+(\d+),\.*'
        pattern2 = r'model:\s+(\d+),\s+ misfit value:\s+(-?\d+\.\d*|-?\d*\.\d+)\.*'
        pattern3 = anywhite+'('+intg+')'+anywhite+real_dn+anywhite+real_dn+anychar
        line = f.readline()

        while len(line) != 0:
            if string.find(line, 'iteration:') != -1:
                match = re.search(pattern1,line)
                iteration = match.group(1)
            if string.find(line, 'model:') != -1:
                self.thickness.append([])
                self.velocity.append([])
                match = re.search(pattern2,line)
                self.misfitval.append(float(match.group(2)))
                for i in range(0,5):
                    newline = f.readline()
                    match = re.search(pattern3,newline)
                    self.velocity[self.ntot].append(float(match.group(3)))
                    self.thickness[self.ntot].append(float(match.group(2)))
                self.ntot = self.ntot + 1

            line = f.readline()
            self.sortedmisfit = sorted(self.misfitval)
     

    def model_transform(self):
        """python version of M. Sambridge's plotting routine 'surface_model_scb.f';
        interpolates between points in model space"""
        for i in range(0,self.ntot):
            
            k=0
            self.depth.append([])
            self.speed.append([])
            self.depth[i].append(0)

            for j in range(0,5):
                try:
                    thick = self.thickness[i][j]
                except Exception, e:
                    print i, self.thickness[i]
                    sys.exit(1)
                if thick > 0.0:
                    v1 = self.velocity[i][j]
                    v2 = v1 + 0.01
                    gr = (v2-v1)/thick

                    nsub =int(thick/2.)+1
                    dh   = thick/float(nsub)

                    for m in range(1,nsub+1):
                        k = k+1
                        self.depth[i].append(self.depth[i][k-1]+dh)
                        self.speed[i].append(v1 + gr*dh*float(m-1))
                        #print "-->depth is: ", self.depth[i][k]
                        #print "-->velocity is: ", self.speed[i][k-1]

            k = k+1
            if self.depth[i][k-1] <= self.dmax:
                self.depth[i].append(self.dmax)
            else:
                self.depth[i].append(self.depth[i][k-1])

            self.speed[i].append(self.speed[i][k-2])
            self.nsublayer.append(k)

            #print "-->depth is: ", self.depth[i][k]
            #print "-->velocity is: ", self.speed[i][k-1]
            #print "-->number of layers: ",self.nsublayer[i]


    def grd2gmt(self,kosu,kosu2,nd,nv):
        """write out matrix values in columns in order to pipe to gmt"""
        for i in range(0,nd):
            for j in range(0,nv):
                if kosu[i][j] <= 0.001 and kosu2[i][j] <= 0.001:
                    print j*self.dx, i*self.dy,'NaN', 'NaN'
                elif kosu[i][j] <= 0.001 and kosu2[i][j] > 0.001:
                    print j*self.dx, i*self.dy,'NaN', '0.5'
                else:
                    print j*self.dx, i*self.dy, kosu[i][j],'0.5'
                    #print j*self.dx, i*self.dy, math.log10(kosu[i][j]),'0.5'



    def model_count(self):
        """python version of M. Sambridge's plotting routine 'surface_model_scb.f';
        counts number of (interpolated) models running over each gridpoint of
        'kosu1' and 'kosu2'"""
        nd = int((self.dmax - self.dmin)/self.dy)+1
        nv = int((self.vmax - self.vmin)/self.dx)+1
        kosu1 = numpy.zeros((nd+1,nv+1))
        kosu2 = numpy.zeros((nd+1,nv+1))

        
        for i in range(0,self.ntot):

            for k in range(1,self.nsublayer[i]+1):
                id1 = int((self.depth[i][k-1]-self.dmin)/self.dy)+1
                if not id1 > nd:
                    id2 = int((self.depth[i][k]-self.dmin)/self.dy)+1
                    if id2 > nd:
                        id2 = nd
                    iv = int((self.speed[i][k-1]-self.vmin)/self.dx)+1
                    for ii in range(id1,id2+1):
                        if self.misfitval[i] < self.sortedmisfit[self.best]:
                            kosu1[ii][iv] = kosu1[ii][iv]+1
                        kosu2[ii][iv] = kosu2[ii][iv]+1
                            

                if k >= 2:
                    iv1 = int((self.speed[i][k-2]-self.vmin)/self.dx)+1
                    iv2 = int((self.speed[i][k-1]-self.vmin)/self.dx)+1
                    ii  = int((self.depth[i][k-1]-self.dmin)/self.dy)+1
                    if not ii > nd:
                        if iv1 < iv2:
                            for iv in range(iv1,iv2+1):
                                if self.misfitval[i] < self.sortedmisfit[self.best]:
                                    kosu1[ii][iv] = kosu1[ii][iv]+1
                                kosu2[ii][iv] = kosu2[ii][iv]+1
                        elif iv1 > iv2:
                            for iv in range(iv2,iv1-1,-1):
                                if self.misfitval[i] < self.sortedmisfit[self.best]:
                                    kosu1[ii][iv] = kosu1[ii][iv]+1
                                kosu2[ii][iv] = kosu2[ii][iv]+1
                        else:
                            if self.misfitval[i] < self.sortedmisfit[self.best]:
                                kosu1[ii][iv] = kosu1[ii][iv]-1
                            kosu2[ii][iv] = kosu2[ii][iv]-1

        self.grd2gmt(kosu1,kosu2,nd,nv)

if __name__ == '__main__':

    try:
        modelname = sys.argv[1]
    except:
        print "Usage: ", sys.argv[0]," model-filename"; sys.exit(1)
        
    new=DensPlotGrd()
    new.get_models(modelname)
    new.model_transform()
    new.model_count()
