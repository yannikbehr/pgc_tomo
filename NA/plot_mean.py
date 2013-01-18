#!/usr/bin/env python
"""read in final model from sambridge na-inversion and print vs-velocity and depth
to stdout\n
$Author$
$Rev:$
$LastChangedDate:$
"""

import string, optparse, sys

class PlotVsModel:

    def __init__(self):
        self.depth = []
        self.vel   = []
        self.coord = [] 
        self.dmax = 100


    def read_model_file(self, filename):
        f = open(filename,'r')
        f.readline()
        f.readline()
        line = f.readline()
        depth = 0
        while len(line) != 0:
            a = string.split(line)
            depth = depth + 0.5
            self.vel.append(float(a[0]))
            self.depth.append(depth)
            line = f.readline()
        
    def create_coord(self):
        self.coord.append((0,self.vel[0]))
        for i in range(0,len(self.depth)):
            if i < len(self.depth)-1:
                self.coord.append((self.depth[i],self.vel[i]))
                self.coord.append((self.depth[i],self.vel[i+1]))
            elif i == len(self.depth)-1:
                self.coord.append((self.depth[i],self.vel[i]))
        # adding semi-infinite layer
        self.coord.append((self.dmax,self.vel[i]))
        for i,j in self.coord:
            print j,i

    
if __name__ == '__main__':
    cmdargs = []
    if len(sys.argv) < 1:
        cmdargs.append("-h")
    else:
        cmdargs = sys.argv[1:]
    p = optparse.OptionParser()
    p.add_option('--file','-f')
    options, arguments = p.parse_args(args=cmdargs)
    test = PlotVsModel()
    test.read_model_file(options.file)
    test.create_coord()
