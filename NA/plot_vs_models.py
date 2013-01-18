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
        contents = f.readlines()
        for i in range(0,len(contents)):
            if string.find(contents[i], 'Final model')!=-1:
                for j in range(i+2,len(contents)):
                    a = string.split(contents[j])
                    b = string.rstrip(a[2],')')
                    self.depth.append(float(b))
                    self.vel.append(float(a[3]))


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
