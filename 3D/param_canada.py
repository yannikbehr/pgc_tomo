#!/usr/bin/env mypython
"""
Create a dinver compatible xml-file from a given model of layer
thicknesses, compressional and shear velocity and density.
"""

import sys, os, os.path
from xml.dom.minidom import parse
from xml.dom.minidom import getDOMImplementation
import xml.dom.minidom
from subprocess import *
import numpy as np
from gmtpy import GMT

def fixed_writexml(self, writer, indent="", addindent="", newl=""):
    """
    Writes well-formatted xml. The default writexml routine creates
    files that can't be read by dinver.
    """
    # indent = current indentation
    # addindent = indentation to add to higher levels
    # newl = newline string
    writer.write(indent+"<" + self.tagName)
    
    attrs = self._get_attributes()
    a_names = attrs.keys()
    a_names.sort()
    
    for a_name in a_names:
        writer.write(" %s=\"" % a_name)
        xml.dom.minidom._write_data(writer, attrs[a_name].value)
        writer.write("\"")
    if self.childNodes:
        if len(self.childNodes) == 1 and self.childNodes[0].nodeType == 3:
            writer.write(">%s</%s>%s" % (self.childNodes[0].data, self.tagName, newl))
            return
        writer.write(">%s"%(newl))
        for node in self.childNodes:
            node.writexml(writer,indent+addindent,addindent,newl)
        writer.write("%s</%s>%s" % (indent,self.tagName,newl))
    else:
        writer.write("/>%s"%(newl))
        # replace minidom's function with ours
        xml.dom.minidom.Element.writexml = fixed_writexml
                                                                                                                                                            
oldwritexml = xml.dom.minidom.Element.writexml
xml.dom.minidom.Element.writexml = fixed_writexml

class DinverParam:
    """
    Class to create xml files that conform to the dinver requirements.
    """
    
    def __init__(self):
        self.npl = 0
        self.nsl = 0
        self.nrl = 0
        self.nul = 0
        impl = getDOMImplementation()
        self.doc = impl.createDocument(None, "Dinver", None)
        self.docel = self.doc.documentElement

    def addtextentry(self,name,value):
        node = self.doc.createElement(name)
        text = self.doc.createTextNode(value)
        node.appendChild(text)
        return node

    def add_header(self):
        plgtag = self.doc.createElement('pluginTag')
        txtplgtag = self.doc.createTextNode('DispersionCurve')
        plgtag.appendChild(txtplgtag)
        plgttl = self.doc.createElement('pluginTitle')
        txtplgttl = self.doc.createTextNode('Surface Wave Inversion')
        plgttl.appendChild(txtplgttl)
        self.top = self.doc.createElement('ParamGroundModel')
        self.docel.appendChild(plgtag)
        self.docel.appendChild(plgttl)
        self.docel.appendChild(self.top)

    def add_footer(self):
        spscript = self.doc.createElement('ParamSpaceScript')
        sptxt = self.doc.createElement('text')
        text = self.doc.createTextNode("""// linear(&quot;TopVs1&quot;,&quot;&gt;&quot;,1.2,&quot;BottomVs0&quot;,1000);
// which means: TopVs1&gt;1.2*BottomVs0+1000""")
        sptxt.appendChild(text)
        spscript.appendChild(sptxt)
        self.top.appendChild(spscript)

    def add_prm_layer(self,prmlayer,mint,maxt,minb,maxb,linkto,dhmin,dhmax,isDepth,lvz=False):
        """
        Add nodes that are common to all entries.
        """
        prmlayer.appendChild(self.addtextentry('shape','Uniform'))
        if not lvz:
            prmlayer.appendChild(self.addtextentry('lastParamCondition','true'))
        else:
            prmlayer.appendChild(self.addtextentry('lastParamCondition','false'))
        prmlayer.appendChild(self.addtextentry('nSubayers',str(5)))
        prmlayer.appendChild(self.addtextentry('topMin',str(mint)))
        prmlayer.appendChild(self.addtextentry('topMax',str(maxt)))
        prmlayer.appendChild(self.addtextentry('bottomMin',str(minb)))
        prmlayer.appendChild(self.addtextentry('bottomMax',str(maxb)))
        prmlayer.appendChild(self.addtextentry('linkedTo',linkto))
        prmlayer.appendChild(self.addtextentry('isDepth',isDepth))
        prmlayer.appendChild(self.addtextentry('dhMin',str(dhmin)))
        prmlayer.appendChild(self.addtextentry('dhMax',str(dhmax)))
        
    def add_p_layer(self,vmint,vmaxt,vminb,vmaxb,linkto,dhmin=1,dhmax=100,isDepth='true',lvz=False):
        """
        Add compressional velocity layer.
        """
        if self.npl < 1:
            self.pmdl = self.doc.createElement('ParamProfile')
            self.pmdl.appendChild(self.addtextentry('type','Param'))
            self.pmdl.appendChild(self.addtextentry('longName','Compression-wave velocity'))
            self.pmdl.appendChild(self.addtextentry('shortName','Vp'))
            self.pmdl.appendChild(self.addtextentry('unit','m/s'))
            self.pmdl.appendChild(self.addtextentry('defaultMinimum',str(200)))
            self.pmdl.appendChild(self.addtextentry('defaultMaximum',str(5000)))
            self.pmdl.appendChild(self.addtextentry('defaultCondition','LessThan'))
            self.top.appendChild(self.pmdl)

        prmlayer = self.doc.createElement('ParamLayer')
        prmlayer.setAttribute('name','Vp%d'%self.npl)
        self.add_prm_layer(prmlayer,vmint,vmaxt,vminb,vmaxb,linkto,dhmin,dhmax,isDepth,lvz=lvz)
        self.pmdl.appendChild(prmlayer)
        self.npl += 1
        
    def add_nu_layer(self,numint,numaxt,numinb,numaxb,linkto,dhmin=1,dhmax=100,isDepth='true',lvz=False):
        """
        Add Poisson's ratio layer.
        """
        if self.nul < 1:
            self.numdl = self.doc.createElement('ParamProfile')
            self.numdl.appendChild(self.addtextentry('type','Condition'))
            self.numdl.appendChild(self.addtextentry('longName',"Poisson's Ratio"))
            self.numdl.appendChild(self.addtextentry('shortName','Nu'))
            self.numdl.appendChild(self.addtextentry('unit',''))
            self.numdl.appendChild(self.addtextentry('defaultMinimum',str(0.2)))
            self.numdl.appendChild(self.addtextentry('defaultMaximum',str(0.5)))
            self.numdl.appendChild(self.addtextentry('defaultCondition','GreaterThan'))
            self.top.appendChild(self.numdl)

        prmlayer = self.doc.createElement('ParamLayer')
        prmlayer.setAttribute('name','Nu%d'%self.nul)
        self.add_prm_layer(prmlayer,numint,numaxt,numinb,numaxb,linkto,dhmin,dhmax,isDepth,lvz=lvz)
        self.numdl.appendChild(prmlayer)
        self.nul += 1
        
    def add_s_layer(self,vmint,vmaxt,vminb,vmaxb,linkto,dhmin=1,dhmax=100,isDepth='true',lvz=False):
        """
        Add shear velocity layer.
        """
        if self.nsl < 1:
            self.smdl = self.doc.createElement('ParamProfile')
            self.smdl.appendChild(self.addtextentry('type','Param'))
            self.smdl.appendChild(self.addtextentry('longName','Shear-wave velocity'))
            self.smdl.appendChild(self.addtextentry('shortName','Vs'))
            self.smdl.appendChild(self.addtextentry('unit','m/s'))
            self.smdl.appendChild(self.addtextentry('defaultMinimum',str(150)))
            self.smdl.appendChild(self.addtextentry('defaultMaximum',str(3500)))
            self.smdl.appendChild(self.addtextentry('defaultCondition','LessThan'))
            self.top.appendChild(self.smdl)

        prmlayer = self.doc.createElement('ParamLayer')
        prmlayer.setAttribute('name','Vs%d'%self.nsl)
        self.add_prm_layer(prmlayer,vmint,vmaxt,vminb,vmaxb,linkto,dhmin,dhmax,isDepth)
        self.smdl.appendChild(prmlayer)
        self.nsl += 1
        
    def add_rho_layer(self,rmint,rmaxt,rminb,rmaxb,linkto,dhmin=1,dhmax=100,isDepth='true',lvz=False):
        """
        Add density layer.
        """
        if self.nrl < 1:
            self.rmdl = self.doc.createElement('ParamProfile')
            self.rmdl.appendChild(self.addtextentry('type','Param'))
            self.rmdl.appendChild(self.addtextentry('longName','Density'))
            self.rmdl.appendChild(self.addtextentry('shortName','Rho'))
            self.rmdl.appendChild(self.addtextentry('unit','kg/m3'))
            self.rmdl.appendChild(self.addtextentry('defaultMinimum',str(2000)))
            self.rmdl.appendChild(self.addtextentry('defaultMaximum',str(2000)))
            self.rmdl.appendChild(self.addtextentry('defaultCondition','LessThan'))
            self.top.appendChild(self.rmdl)

        prmlayer = self.doc.createElement('ParamLayer')
        prmlayer.setAttribute('name','Rho%d'%self.nrl)
        self.add_prm_layer(prmlayer,rmint,rmaxt,rminb,rmaxb,linkto,dhmin,dhmax,isDepth,lvz=lvz)
        self.rmdl.appendChild(prmlayer)
        self.nrl += 1

        
    def savexml(self,fout):
        """
        Write xml file to disk.
        """
        f = open(fout,'w')
        self.doc.writexml(f,encoding='UTF-8',addindent='  ',newl='\n')
        f.close()

def plotmdl(th,vs,crth,dplith,lat,lon,prem=False):
    """
    Plot the model constructed from crust2.0 and appended by PREM.
    """
    d = 0
    depth = []
    vels = []
    for t,v in zip(th,vs):
        depth.append(-d/1000.)
        vels.append(v/1000.)
        d += t
        depth.append(-d/1000.)
        vels.append(v/1000.)
    plot(vels,depth,'k',label='Combined model')
    hlines(-crth,1.,6.,'b',label='Crust 2.0')
    hlines(-dplith,1.,6.,'r',label='Lith 5.0')
    ymin,ymax = ylim()
    
    if prem:
        th = [12.0,10.0,15.,20.,20.,35.,35.,35.,35.,45.,45.,60.,30.,50.,50.,50.,50.,35.,35]
        vs = [3.191,3.889,4.479,4.473,4.465,4.457,4.363,4.350,4.338,
              4.325,4.620,4.651,4.683,4.714,5.019,5.163,5.307,5.451,5.478]
        d = 0
        depth = []
        vels = []
        for t,v in zip(th,vs):
            depth.append(-d)
            vels.append(v)
            d += t
            depth.append(-d)
            vels.append(v)
        plot(vels,depth,'k--',label='PREM')
    title('Lat: %.2f  Lon: %.2f'%(lat,lon))
    legend(loc='upper left')
    xlim(1.,6.)
    ylim(ymin,ymax)
    ylabel('Depth [km]')
    xlabel('Vs [km/s]')
    
        
def mdl2param(th,vp,vs,rho,fout):
    """
    Create a paramter file that can be read by dinver.
    """
    param = DinverParam()
    param.add_header()
    cnt = 0
    for _th,_vp,_vs,_rho in zip(th[0:-1],vp[0:-1],vs[0:-1],rho[0:-1]):
        fct = 0.2
        vsmin = _vs - fct * _vs
        vsmax = _vs + fct * _vs
        vpmin = _vp - fct * _vp
        vpmax = _vp + fct * _vp
        thmin = _th - fct * _th
        thmax = _th + fct * _th
        if vsmin <= 0.: vsmin = _vs
        if vpmin <= 0.: vpmin = _vp
        if thmin <= 0.: thmin = _th
        
        param.add_p_layer(vpmin,vpmax,vpmin,vpmax,'Vs%d'%cnt)
        param.add_nu_layer(0.2,0.5,0.2,0.5,'Vs%d'%cnt)
        param.add_s_layer(vsmin,vsmax,vsmin,vsmax,'Not linked',dhmin=thmin,dhmax=thmax,isDepth='false')
        param.add_rho_layer(_rho,_rho,_rho,_rho,'Vs%d'%cnt,lvz=True)
        cnt += 1
    _th, _vp, _vs, _rho = (th[-1],vp[-1],vs[-1],rho[-1])
    vsmin = _vs - fct * _vs
    vsmax = _vs + fct * _vs
    vpmin = _vp - fct * _vp
    vpmax = _vp + fct * _vp
    param.add_p_layer(vpmin,vpmax,vpmin,vpmax,'Not linked')
    param.add_nu_layer(0.2,0.5,0.2,0.5,'Not linked')
    param.add_s_layer(vsmin,vsmax,vsmin,vsmax,'Not linked')
    param.add_rho_layer(_rho,_rho,_rho,_rho,'Not linked',lvz=True)
    param.add_footer()
    param.savexml(fout)


def addprem(thcrust):
    """
    Extend the crustal model from crust2.0 by 100 km with values from
    PREM. The crustal layers of PREM get ignored by setting the
    velocities of the crustal layers to those of the uppermost mantle
    layer. This avoids artificial low velocity zones beneath the Moho
    because the crustal thickness of Crust2.0 is smaller than that of
    PREM.
    """
    th   = [12.0,10.0,15.,20.,20.,35.,35.,35.,35.,45.,45.,60.,30.,50.,50.,50.,50.,35.,35]
    ### Original PREM values
    #rho = [2.6,2.9,3.381,3.379,3.377,3.375,3.371,3.367,3.363,3.436,
    #       3.463,3.49,3.516,3.543,3.787,3.850,3.913,3.976,3.984]
    #vp  = [5.793,6.792,8.101,8.091,8.079,8.067,7.984,7.963,7.942,
    #       7.920,8.606,8.692,8.778,8.865,9.347,9.601,9.856,10.111,10.165]
    #vs  = [3.191,3.889,4.479,4.473,4.465,4.457,4.363,4.350,4.338,
    #       4.325,4.620,4.651,4.683,4.714,5.019,5.163,5.307,5.451,5.478]
    ### PREM values with the crustal layers set to upper mantle velocities and density
    rho = [3.381,3.381,3.381,3.379,3.377,3.375,3.371,3.367,3.363,3.436,
           3.463,3.49,3.516,3.543,3.787,3.850,3.913,3.976,3.984]
    vp  = [8.101,8.101,8.101,8.091,8.079,8.067,7.984,7.963,7.942,
           7.920,8.606,8.692,8.778,8.865,9.347,9.601,9.856,10.111,10.165]
    vs  = [4.479,4.479,4.479,4.473,4.465,4.457,4.363,4.350,4.338,
           4.325,4.620,4.651,4.683,4.714,5.019,5.163,5.307,5.451,5.478]
    d = 0
    dmin = 1e10
    for i,h in enumerate(th):
        d += h
        if d > thcrust and (d-h) < thcrust:
            imin = i
            dmin = d
        elif d == thcrust:
            imin = i + 1
            dmin = d + th[i+1]
        if thcrust+100 < d and thcrust+100 > d-h:
            imax = i
            dmax = d
            break
    th[imin] = dmin-thcrust
    th[imax] = thcrust+100-(dmax-h)
    
    return np.array(th[imin:imax+1])*1000, np.array(vp[imin:imax+1])*1000,\
           np.array(vs[imin:imax+1])*1000, np.array(rho[imin:imax+1])*1000

def lith5mdl(lat,lon):
    """
    Get moho depth from Lith5.0 model.
    """
    lons,lats,dp = np.loadtxt('./lith5.0_moho.txt',unpack=True)
    try:
        ind = np.where((lons == lon) & (lats == lat))[0][0]
    except IndexError:
        print "point at lat=%f, lon=%f is out of range of Lith5.0"%(lat,lon)
        print "Crust2.0 is used instead."
        return None
    else:
        if dp[ind] > 1000. or np.isnan(dp[ind]):
            return None
        else:
            return dp[ind]

def crust2mdl(lat,lon,fout,lith50=False):
    """
    Extract the model for a given latitude and longitude from
    Crust2.0. As an alternative to the crustal thickness from Crust2.0
    the moho depth from the Lith5.0 can be used.
    """
    dplith = lith5mdl(lat,lon)
    if dplith is None:
        lith50 = False
    else:
        mohodp = dplith
    rootdir='./crust2'
    curdir = os.getcwd()
    os.chdir(rootdir)
    fcrust = 'outcr'
    if os.path.isfile(fcrust):
        os.remove(fcrust)
    p = Popen(['./getCN2point >/dev/null'],shell=True,stdin=PIPE)
    print >>p.stdin, '%f, %f'%(lat,lon)
    print >>p.stdin, '*'
    p.stdin.close()
    sts = os.waitpid(p.pid, 0)[1]
    f = open(fcrust)
    while True:
        line = f.readline()
        if not line: break
        if line.find('crustal thickness, ave. vp, vs, rho:') != -1:
            a = line.split()
            crth = float(a[6])
        if line.find('layer crustal model') != -1:
            th = []
            vp = []
            vs = []
            rho = []
            d = 0.
            if lith50:
                ### Moho depth is taken from Lith5.0
                while True:
                    line = f.readline()
                    if not line: break
                    t,p,s,r = map(float,line.split()[0:4])
                    if t < 0.01: continue
                    if s < 0.01: continue
                    d += t
                    if d < dplith:
                        th.append(t*1000)
                        vp.append(p*1000)
                        vs.append(s*1000)
                        rho.append(r*1000)
                    else:
                        t -= d-dplith
                        if t > 0.:
                            th.append(t*1000)
                            vp.append(p*1000)
                            vs.append(s*1000)
                            rho.append(r*1000)
            else:
                mohodp = crth
                while True:
                    line = f.readline()
                    if not line: break
                    t,p,s,r = map(float,line.split()[0:4])
                    if t < 0.01: continue
                    if s < 0.01: continue
                    th.append(t*1000)
                    vp.append(p*1000)
                    vs.append(s*1000)
                    rho.append(r*1000)
                        
    os.chdir(curdir)
    thh, vpp, vss, rhoo = addprem(mohodp)
    if 0:
        ### enable for debugging only
        print crth
        for t,p,s,r in zip(th,vp,vs,rho):
            print t, p, s, r
        if dplith is None:
            dplith = 0.
        plotmdl(np.r_[np.array(th),thh],np.r_[np.array(vs),vss],crth,dplith,lat,lon,prem=True)
        
    mdl2param(np.r_[np.array(th),thh],np.r_[np.array(vp),vpp],np.r_[np.array(vs),vss],np.r_[np.array(rho),rhoo],fout)
    #mdl2param(np.array(th),np.array(vp),np.array(vs),np.array(rho),fout)
    


if __name__ == '__main__':
    from pylab import *
    fout = 'canada_45_-75.param'
    if 0:
        runlat = arange(24.,89.)
        runlon = arange(-156.,-43.)
        for lat in runlat:
            for lon in runlon:
                crust2mdl(lat,lon,fout)
    if 1:
        crust2mdl(35,-52,fout,lith50=True)
                
    #dplith = lith5mdl(-45,-75)
    #print dplith
    
