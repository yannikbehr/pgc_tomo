#!/usr/local/bin/python


"""run NA inversion on evison-grid at GNS"""


import sys, os, shutil
from os.path import basename,join,isfile
sys.path.append('./py_ext/na/writenad')
sys.path.append('./py_ext/na/density')
sys.path.append(join(os.environ['PROC_SRC'],'disp_curves'))
import mywritenad as wnad
import writenad as w
from pylab import *
import convdisp
from nadplot import *
from gmtpy import GMT
from scipy import optimize

def rmkdir(dirn):
    cmd = "ssh -i %s %s 'mkdir %s' "%(key,host,dirn)
    print cmd;err=os.system(cmd)
    if err != 0 and err != 256: return 0
    return 1


def movedata(fn,dest):
    cmd = "scp -i %s %s %s:%s/%s"%(key,fn,host,dest,basename(fn))
    print cmd;err=os.system(cmd)
    if err != 0: return 0
    return 1


def runscript(name,rdir):
    cmd = "ssh -i %s %s 'cd %s && tcsh %s&' "%(key,host,dest,basename(rdir)+'/'+basename(name))
    print cmd;err=os.system(cmd)
    if err != 0: return 0
    return 1


def getres(fn,dest):
    cmd = "scp -i %s %s:%s %s/%s"%(key,host,fn,dest,basename(fn))
    print cmd;err=os.system(cmd)
    if err != 0: return 0
    return 1
    

def writescript(scname,dispfn,dsptype,wtype,rno,velmdl,st1,st2):
    f=open(scname,'w')
    myvars = {}
    myvars['obsdisp']  = dispfn
    myvars['dsptype']  = dsptype
    myvars['wtype']    = wtype
    myvars['runnum']   = rno
    myvars['velmodel'] = velmdl
    myvars['stat1']    = st1
    myvars['stat2']    = st2
    f.write("""
set inv = %(wtype)s 
if ( $inv == 'ray' ) then
  set dispcurve_rayleigh = %(obsdisp)s
endif
if ( $inv == 'love' ) then
  set dispcurve_love = %(obsdisp)s
endif
set WaveType = %(dsptype)s
if ( $WaveType == "phase" ) set porg = 1
if ( $WaveType == "group" ) set porg = 2
set runnum = %(runnum)d
set velmodel = %(velmodel)s
set station1 = %(stat1)s
set station2 = %(stat2)s

if ( ! -e ./summary.$runnum ) echo " " > ./summary.$runnum 
 
echo " " >> summary.$runnum 
set date1 = `date | awk '{print $1, $2, $3, $4}'`
echo This ran on :  $date1 >> summary.$runnum
echo " " >> summary.$runnum 
echo "Input velocity model is $velmodel" 
echo "Run number is $runnum" 
echo " " >> summary.$runnum 
echo " " >> summary.$runnum 
echo Velocity-model $velmodel >> summary.$runnum
echo " " >> summary.$runnum 
more ./models/$velmodel >> summary.$runnum
more na.in >> summary.$runnum
echo " " >> summary.$runnum 
echo "observed dispersion curve:" >> summary.$runnum

switch (${inv})

case "joint":
    echo "Joint inversion of Rayleigh and Love waves"
    if ( -e surface.in ) rm  surface.in
    echo $velmodel | awk -f form_surface_in.awk > surface.in

    set LorRorJ = 3
    echo $LorRorJ >> surface.in
    echo $porg >> surface.in

    set nwaves = 2
    echo $nwaves >> surface.in
    rm sobs.d
    cp sobs_master.d sobs.d
    echo $dispcurve_rayleigh >> sobs.d
    echo $dispcurve_love >> sobs.d
    echo $dispcurve_rayleigh >> surface.in
    echo $dispcurve_love >> surface.in
    echo $dispcurve_rayleigh >> summary.$runnum
    echo $dispcurve_love >> summary.$runnum

breaksw

case "ray":
    echo "Rayleigh wave inversion"
    set LorRorJ = 2
    rm  surface.in
    echo $velmodel | awk -f form_surface_in.awk > surface.in
    echo $LorRorJ >> surface.in
    echo $porg >> surface.in

    set nwaves = 1
    echo $nwaves >> surface.in
    rm sobs.d
    cp sobs_master.d sobs.d
    echo ${dispcurve_rayleigh} >> sobs.d
    echo ${dispcurve_rayleigh} >> surface.in

    
    echo $dispcurve_rayleigh >> summary.$runnum

breaksw

case "love":
    echo "Love wave inversion"
    set LorRorJ = 1
    rm  surface.in
    echo $velmodel | awk -f form_surface_in.awk > surface.in
    echo $LorRorJ >> surface.in
    echo $porg >> surface.in
    set nwaves = 1
    echo $nwaves >> surface.in
    rm sobs.d
    cp sobs_master.d sobs.d
    echo ${dispcurve_love} >> sobs.d
    echo ${dispcurve_love} >> surface.in
    echo $dispcurve_love >> summary.$runnum
breaksw

default:
    echo "Joint inversion of Rayleigh and Love waves"
breaksw

endsw

more surface.in

set modeldir = ./models
if ( -e ${modeldir}/models ) rm ${modeldir}/models

echo "starting surface_na"
#../bin/surface_na
../bin/surface_na_geopsy

# finish off the summary file
# ===========================

mv summary.$runnum ./sumfiles/summary.$runnum.$station1.$station2
if ( -e na.nad ) mv na.nad ./nadfiles/na-${velmodel}-r${runnum}.nad

#cp fort.50  ${velmodel}-r${runnum}-calculated.dat
#cp fort.51  observed-trace1.dat

if (  ! -e ${modeldir}/models ) then
    echo "no models file"
    exit
endif 

set modeldir = ./models
cp ${modeldir}/${velmodel}   ${modeldir}/surface_param-r${runnum}
cp ${modeldir}/models  ${modeldir}/models-${velmodel}-${station1}-${station2}-r${runnum}
mv fort.60 ./${modeldir}/fort-${velmodel}-${station1}-${station2}-r${runnum}.60
mv fort.61 ./${modeldir}/fort-${velmodel}-${station1}-${station2}-r${runnum}.61
"""%myvars)
    f.close()
    return 1


def writenain(fnain,iter=100,sinit=50,sample=100,ncells=50):
    f = open(fnain,'w')
    f.write("""#
# Neighbourhood Algorithm input options file
#
0           : Algorithm type (NA or Uniform MC: 1=MC,0=NA)
%(iter)d    : Maximum number of iterations
%(sinit)d   : Sample size for first iteration
%(sample)d  : Sample size for all other iterations
%(ncells)d  : Number of cells to re-sample
n,210728    : Use Quasi random number generator ? (y/n);random seed
0           : Type of initial sample (0=random;1=read in a NAD file)
1           : Output information level (0=silent,1=summary info,2=1+models)
y           : Turn timing mode on ? (y/n)
n           : Turn debug mode on ? (y/n)
"""%vars())
    return 1
            
    

def txt2nad(fn1,fn2,fnad):
    models,data,ns1,ns,nr,nd,ne = wnad.get_models(fn1)
    ranges, scales = wnad.get_ranges(fn2,nd)
    newmodel = zeros((nd,ne+1))
    nhmax = 10000
    header = zeros(nhmax)
    for jj in range(ne+1):
        for kk in range(nd):
            newmodel[kk,jj]=models[jj*nd+kk]

    lu = 16
    nh = 256
    fout = fnad
    tmpfile = 'junk'
    w.write_header(lu,tmpfile,nh,ranges,scales,ns1,ns,nr)
    header = w.read_header(lu,tmpfile,nh,nd)
    nhu = 0
    iform = 1
    w.write_nad(lu,fout,nhu,header[0:nh],iform,newmodel,data[0:ne+1],nd,ne+1,nh)
    if isfile(tmpfile):
        os.remove(tmpfile)
    return 1


def plotnad(fnad,fout):
    gmt = GMT(config={'BASEMAP_TYPE':'plain','ANOT_FONT_SIZE':8,
                      'LABEL_FONT_SIZE':10,'COLOR_BACKGROUND':'255/255/255',
                      'COLOR_FOREGROUND':'0/0/0','COLOR_NAN':'255/255/255'} )

    grd=gmt.tempfilename('tmp.grd')
    grdcpt=gmt.tempfilename('tmp.cpt')
    xyz=gmt.tempfilename('xyz.txt')
    xyz2=gmt.tempfilename('xyz2.txt')
    rng='1/5/0/40'
    scl='X4.2/-6'
    kosu1,kosu2,x,y,dbest,sbest,dmean,smean = nadplot(fnad,smin=1.5)
    f = open(xyz,'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if kosu1[ii,jj] > 0:
                print >>f, x[jj], y[ii], kosu1[ii,jj]
    f.close()
    f = open(xyz2,'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if kosu2[ii,jj] > 0:
                print >>f, x[jj], y[ii], '0.5'

    f.close()
    gmt.xyz2grd(xyz,G=grd,R=rng,I='.02/.5',out_discard=True)
    gmt.grd2cpt(grd,C="wysiwyg",Q='o',Z=True,out_filename=grdcpt)
    gmt.psmask(xyz2,R=rng,T=True,J=scl,I='.02/.5',G='lightgray')
    gmt.grdimage(grd,J=scl,R=rng,Q=True,P=True,C=grdcpt)
    gmt.psbasemap(R=rng,J=scl,B='a1f.5:Velocity [km/s]:/a10f5:Depth [km]:WNse')
    gmt.psxy(R=True,J=True,B=True,W='3,red',in_columns=[sbest,dbest])
    gmt.psxy(R=True,J=True,B=True,W='3,white',in_columns=[smean,dmean])
    gmt.psscale(C=grdcpt,D='1.0/1./4c/.4ch',B='a40f10:No. of models:/::')
    gmt.save(fout) 
    os.system('gv '+fout+'&')
    return 1


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
                      'COLOR_FOREGROUND':'0/0/0','COLOR_NAN':'255/255/255'} )
    data,ok = read_nab(nabfile)
    if not ok: return 0
    xshift = 1
    yshift = 2
    for ii in range(0,4):
        x = data[ii,0,:]; y = data[ii,1,:]
        rng = '%f/%f/%f/%f'%(floor(x[0]),ceil(x[-1]),y.min(),y.max()+0.1*y.max())
        scl = 'X4c/4c'
        ant = (ceil(x[-1])-floor(x[0]))/5.
        fnt = ant/2
        gmt.psxy(R=rng,J=scl,B='a%ff%f:Layer thickness [km]:weSn'%(ant,fnt),W='3,black',X='a%dc'%xshift,Y='a%dc'%yshift,in_columns=[x,y])
        xshift = xshift + 5
    yshift = 8.5
    xshift = 1
    for ii in range(5,9):
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


    
def run_bayes(nadfile,nabfile,nabin,nabplotf):
    if not write_nabin(nabin,nstep=200): return 0
    if not os.path.isfile(nabin): return 0
    baydir = '/home/behrya/dev/na_bayes_svn/Demo/data/'
    shutil.copy(nabin,baydir+'nab.in')
    shutil.copy(nadfile,baydir+'ensemble.nad')
    os.chdir(baydir)
    os.system('./surf_nab')
    shutil.copy('nab.out',nabfile)
    if not ok: return 0
    return 1


if __name__ == '__main__':
    ##### parameters #####################################
    #ftanfile = '/data/sabine/yannik/Results/stack/horizontal/nord/COR_WCZ_TIKO.SAC_TT_s_2_DISP.3'
    ftanfile = '/data/sabine/yannik/Results/stack/horizontal/nord/COR_MATA_TIKO.SAC_TT_s_2_DISP.4'
    filename = '/home/behrya/dev/data/tiko_mata_na_l_tt.txt'
    scname   = '../data/na_remote.csh'
    fnain    = '../data/na.in'
    key      = '/home/behrya/.ssh/evison_dsa'
    host     = 'behrya@evison.gns.cri.nz'
    dest     = '/home/behrya/sambridge_with_love_w/data/'
    #rdir     = '/home/behrya/sambridge_with_love_w/data/run_remote'
    rdir     = '/home/behrya/sambridge-geopsy_svn/data/run_remote'
    resdir   = '/home/behrya/dev/data/'
    dispfn   = join(basename(rdir),basename(filename))
    wtype    = 'love'
    dsptype  = 'group'
    rno      = 7006
    velmdl   = 'northland13'
    st1      = 'tiko'
    st2      = 'mata'
    mdlfn    = join(dest,'models/models-'+velmdl+'-'+st1+'-'+st2+'-r'+str(rno))
    prmfn    = join(dest,'models/surface_param-r'+str(rno))
    mdldsp   = join(dest,'models/fort-'+velmdl+'-'+st1+'-'+st2+'-r'+str(rno)+'.60')
    nadfile  = join(resdir,'models-'+velmdl+'-'+st1+'-'+st2+'-r'+str(rno)+'.nad')
    nabfile  = join(resdir,'models-'+velmdl+'-'+st1+'-'+st2+'-r'+str(rno)+'.nab')
    pdfout   = join(resdir,'models-'+velmdl+'-'+st1+'-'+st2+'-r'+str(rno)+'.pdf')
    nabin    = join(resdir,'nabin_r'+str(rno))
    nabplotf = join(resdir,'nab_1d_r'+str(rno))+'.pdf'
    nabpdfmoho = join(resdir,'nab_1d_moho_r'+str(rno))+'.pdf'
    ######################################################

    #if not writenain(fnain): sys.exit(1)
    #if not convdisp.ConvDisp(ftanfile,filename,wtype,dsptype)(): sys.exit(1)
    if not rmkdir(rdir): sys.exit(1)
    #if not movedata(filename,rdir): sys.exit(1)
    #if not movedata(fnain,dest): sys.exit(1)
    #if not writescript(scname,dispfn,dsptype,wtype,rno,velmdl,st1,st2): sys.exit(1)
    if not movedata(scname,rdir): sys.exit(1)
    if not runscript(scname,rdir): sys.exit(1)
    #if not getres(mdlfn,resdir): sys.exit(1)
    #if not getres(prmfn,resdir): sys.exit(1)
    #if not getres(mdldsp,resdir): sys.exit(1)
    #f1 = join(resdir,basename(mdlfn))
    #f2 = join(resdir,basename(prmfn))
    #if not txt2nad(f1,f2,nadfile): sys.exit(1)
    #if not plotnad(nadfile,pdfout): sys.exit(1)
    #if not run_bayes(nadfile,nabfile,nabin,nabplotf): sys.exit(1)
    #if not plot_moho(nabfile,nabpdfmoho): sys.exit(1)
    #if not plot_pdf(nabfile,nabplotf): sys.exit(1)
