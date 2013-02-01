#!/usr/bin/env python
"""
extract xml-contents from dinver *.param files and write it
as an ascii table to stdout
"""

import sys
import os
import shutil
from xml.dom.minidom import parse

class ProFile:
    def __init__(self, fn=None, factor=1000.0):
        self.vsmin = []
        self.vsmax = []
        self.vpmin = []
        self.vpmax = []
        self.thmin = []
        self.thmax = []
        self.numin = []
        self.numax = []
        self.rhomin = []
        self.rhomax = []
        self.factor = factor
        if fn:
            self.__call__(fn)

    def __call__(self, fname):
        try:
            if os.system("tar -xzf %s" % fname):
                raise Exception('system call failed')
        except Exception, e:
            print "cannot un-tar file!"
            print e
            shutil.copy(fname, 'contents.xml')
        # parse the document
        f = open('contents.xml')
        dom = parse(f)

        for _n in dom.getElementsByTagName('ParamProfile'):
            if _n.getElementsByTagName('shortName')[0].firstChild.data == 'Vs':
                self.vs(_n)
            if _n.getElementsByTagName('shortName')[0].firstChild.data == 'Vp':
                self.vp(_n)
            if _n.getElementsByTagName('shortName')[0].firstChild.data == 'Rho':
                self.rho(_n)
            if _n.getElementsByTagName('shortName')[0].firstChild.data == 'Nu':
                self.nu(_n)


    def __del__(self):
        # ## cleaning up
        if os.path.isfile('contents.xml'):
            os.remove('contents.xml')


    def vs(self, nd):
        layers = nd.getElementsByTagName('ParamLayer')
        for _l in layers:
            self.vsmin.append(float(_l.getElementsByTagName('topMin')[0].firstChild.data) / self.factor)
            self.vsmax.append(float(_l.getElementsByTagName('topMax')[0].firstChild.data) / self.factor)
            self.thmin.append(float(_l.getElementsByTagName('dhMin')[0].firstChild.data) / self.factor)
            self.thmax.append(float(_l.getElementsByTagName('dhMax')[0].firstChild.data) / self.factor)
        self.thmin[-1] = 0.0
        self.thmax[-1] = 0.0
        self.nlayers = len(self.thmin)


    def vp(self, nd):
        layers = nd.getElementsByTagName('ParamLayer')
        for _l in layers:
            self.vpmin.append(float(_l.getElementsByTagName('topMin')[0].firstChild.data) / self.factor)
            self.vpmax.append(float(_l.getElementsByTagName('topMax')[0].firstChild.data) / self.factor)


    def nu(self, nd):
        layers = nd.getElementsByTagName('ParamLayer')
        for _l in layers:
            self.numin.append(float(_l.getElementsByTagName('topMin')[0].firstChild.data))
            self.numax.append(float(_l.getElementsByTagName('topMax')[0].firstChild.data))


    def rho(self, nd):
        layers = nd.getElementsByTagName('ParamLayer')
        for _l in layers:
            self.rhomin.append(float(_l.getElementsByTagName('topMin')[0].firstChild.data) / self.factor)
            self.rhomax.append(float(_l.getElementsByTagName('topMax')[0].firstChild.data) / self.factor)


if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except:
        print "usage %s filename" % os.path.basename(sys.argv[0])
        sys.exit(0)

    prf = ProFile(fname)

    #### print ascii table to sdout
    print "# number of layers: %d" % (prf.nlayers)
    print "#%-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s"\
          % ('Thmin', 'Thmax', 'Vsmin', 'Vsmax', 'Vpmin', 'Vpmax', \
            'Rhomin', 'Rhomax', 'Numin', 'Numax')
    for i in range(len(prf.thmin)):
        print "%-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f"\
              % (prf.thmin[i], prf.thmax[i], prf.vsmin[i], prf.vsmax[i], \
                prf.vpmin[i], prf.vpmax[i], prf.rhomin[i], prf.rhomax[i], \
                prf.numin[i], prf.numax[i])

