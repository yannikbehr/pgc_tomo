#!/usr/local/bin/python
"""
run sambridge's bayesian evaluator on dinver output
"""
from gpdc2nad import *
import shutil, os, sys, os.path
from gmtpy import GMT
from scipy import optimize


def write_nabin(nabin,nrnd=10,nstep=20,marg1d='y',marg2d='n',
                cov='n',fout=50):
    f = open(nabin,'w')
    f.write("""#
#	Input options for program NAB
#
0      :algorithm type(0=NA-resampling,1=uniform Monte Carlo integration)
%(nrnd)d      :number of random walks to perform (each from a different starting cell)
%(nstep)d     :number of steps per random walk 
1      :starting cell of first random walk (1=best fit model)
Q      :use quasi or pseudo-random deviates ? (Q=quasi,P=pseudo)
230728 :seed for pseudo-random deviates, or coefficients in quasi data table.
1      :number of steps before accepting model (in random walk)
2000   :number of steps before refreshing distance list (avoids error build up)
n      :set cpu timing mode
10     :frequency at which diagnostic message is written to standard out 
1      :I/O mode (1=silent,2=verbose,3=debug)
        
        Options for numerical integration

%(marg1d)s      :calculate 1D marginal pdfs ? (Y=yes,n=no)
%(marg2d)s      :calculate 2D marginal pdfs ? (Y=yes,n=no)
%(cov)s      :calculate covariance matrix ? (Y=yes,n=no)
40     :number of bins per axis in 1D marginal pdfs.
20     :number of bins per axis in 2D marginal pdfs.
1     :number of user supplied special functions to integrate (0=none)
0     :number of 2D marginals to calculate, (n2d)
32,13  : first pair of variables for 2D marginal (n2d pairs expected)
33,14  : second pair of variables for 2D marginal
34,15  : plot 3
35,16  : plot 4
36,17  : plot 5
1,7    : plot 6
1,13   : plot 7
2,8    : plot 8
2,14   : plot 9
8,13   : plot 10
19,20  : plot 11
20,21  : plot 12
21,22  : plot 13
22,23  : plot 14
26,12  : plot 15
7,13   : plot 16
8,14   : plot 17
9,15   : plot 18
10,16  : plot 19
11,17  : plot 20
12,18  : plot 21
26,25  : plot 22
%(fout)d :frequency of output to file (in number of samples)
"""%vars())
    return 1

def read_nab(nabfile):
    f = open(nabfile,'r')
    lines = f.readlines()
    ir = 0
    for l in lines:
        if l.find('Number of dim') != -1:
            ndim = int(l.split()[4])
        if l.find('Number of user') != -1:
            nsp = int(l.split()[7])
        if l.find('Number of bins per axis for 1D') != -1:
            ndis = int(l.split()[9])
        ### find last 'Results' entry
        if l.find('Results') != -1:
            ir = lines.index(l,ir+1)
    data = zeros((ndim+nsp,2,ndis))
    ### reading results
    il = 0; ic = 0
    for ii in range(ir,len(lines),1):
        if lines[ii].find('Marginal') != -1:
            for ent in lines[ii+2:ii+2+ndis]:
                data[il,0,ic],data[il,1,ic]=map(float,ent.split())
                ic = ic + 1
            il = il + 1
            ic = 0
    return data,1


def plot_moho(nabfile,fout,gauss=True):
    data,ok = read_nab(nabfile)
    if not ok: return 0
    x = data[20,0,:]; y = data[20,1,:]
    gmt = GMT(config={'BASEMAP_TYPE':'plain','ANOT_FONT_SIZE':8,
                      'LABEL_FONT_SIZE':10,'COLOR_BACKGROUND':'255/255/255',
                      'COLOR_FOREGROUND':'0/0/0','COLOR_NAN':'255/255/255'} )
    rng = '%f/%f/%f/%f'%(floor(x[0]),ceil(x[-1]),y.min(),y.max()+0.1*y.max())
    scl = 'X6c/4c'
    ant = int((ceil(x[-1])-floor(x[0]))/10.)
    fnt = ant/2
    gmt.psxy(R=rng,J=scl,B='a%ff%fweSn'%(ant,fnt),W='3,red',in_columns=[x,y])
    if gauss:
        ### fit a gauss-curve to the pdf
        fitfunc = lambda p,x: (1/sqrt(2*pi*p[0]**2))*exp(-(x-p[1])**2/(2*p[0]**2))
        errfunc = lambda p,x,y: fitfunc(p,x)-y
        gaussian = lambda m,s,x: (1/sqrt(2*pi*s**2))*exp(-(x-m)**2/(2*s**2))
        p0 = [10.,50.]
        p1, success = optimize.leastsq(errfunc,p0[:],args=(x,y))
        sigm,mean = p1
        z = gaussian(mean,sigm,x)
        gmt.psxy(R=rng,J=scl,B='a%ff%fweSn'%(ant,fnt),W='3,black,- -',in_columns=[x,z])
    gmt.save(fout)
    os.system('gv '+fout+'&')
    return 1


def plot_pdf(nabfile,fout,gauss=True):
    gmt = GMT(config={'BASEMAP_TYPE':'plain','ANOT_FONT_SIZE':8,
                      'LABEL_FONT_SIZE':10,'COLOR_BACKGROUND':'255/255/255',
                      'COLOR_FOREGROUND':'0/0/0','COLOR_NAN':'255/255/255',
                      'PAGE_ORIENTATION':'landscape'} )
    data,ok = read_nab(nabfile)
    if not ok: return 0
    xshift = 1
    yshift = 2
 #   for ii in range(0,4):
 #       x = data[ii,0,:]; y = data[ii,1,:]
 #       ### skip layers with zero thickness:
 #       if x[0] == x[-1]:continue
 #       rng = '%f/%f/%f/%f'%(floor(x[0]),ceil(x[-1]),y.min(),y.max()+0.1*y.max())
 #       scl = 'X4c/4c'
 #       ant = (ceil(x[-1])-floor(x[0]))/5.
 #       fnt = ant/2
 #       gmt.psxy(R=rng,J=scl,B='a%ff%f:Layer thickness [km]:weSn'%(ant,fnt),W='3,black',X='a%dc'%xshift,Y='a%dc'%yshift,in_columns=[x,y])
 #       xshift = xshift + 5
 #   yshift = 8.5
    xshift = 1
    for ii in range(10,15):
        x = data[ii,0,:]; y = data[ii,1,:]
        rng = '%f/%f/%f/%f'%(floor(x[0]),ceil(x[-1]),y.min(),y.max()+0.1*y.max())
        scl = 'X4c/4c'
        ant = (ceil(x[-1])-floor(x[0]))/5.
        fnt = ant/2
        gmt.psxy(R=rng,J=scl,B='a%ff%f:S-velocity [km/s]:weSn'%(ant,fnt),W='3,black',X='a%dc'%xshift,Y='a%dc'%yshift,in_columns=[x,y])
        xshift = xshift + 5
    gmt.save(fout)
    os.system('gv '+fout+'&')
    return 1


    
def run_bayes(nadfile,nabfile,nabin):
    write_nabin(nabin,nrnd=600,nstep=100,cov='y')
    os.path.isfile(nabin)
    baydir = '/home/behrya/dev/na_bayes_svn/Demo/data/'
    shutil.copy(nabin,baydir+'nab.in')
    shutil.copy(nadfile,baydir+'ensemble.nad')
    os.chdir(baydir)
    os.system('./surf_nab')
    shutil.copy('nab.out',nabfile)



def main(mdlname):
    nadf = mdlname.split('.report')[0]+'.nad'
    nabf = mdlname.split('.report')[0]+'.nab'
    nabin = mdlname.split('.report')[0]+'.nabin'
    ns=100; ns1=50; nit=600
    gpdc2nad(mdlname,prmf,nadf,ns,ns1,nit)
    run_bayes(nadf,nabf,nabin)
    plot_pdf(nabf,'test.ps',gauss=False)
    plot_moho(nabf,'test_moho.ps',gauss=False)


if __name__ == '__main__':
    mdlname = '/home/behrya/dev/data/mt_fixed_moho_ray_c_u_reports/mt_fixed_moho_ray_c_u_01.report'
    prmf = '/home/behrya/dev/data/models/fixed_moho.param'
    mdlname = '/home/behrya/dev/data/mt_ray_u_c_reports/mt_ray_u_c_04.report'
    prmf = '/home/behrya/dev/data/models/northland13_3.param'
    main(mdlname)
