#!/usr/bin/env mypython
"""
Plot results of 1D inversion at discrete geographical points.
"""

# from pylab import *
from scipy import optimize
# from matplotlib import rcParams
import tempfile
import sys
import os
sys.path.append('../NA')
from gpdcreport2mdl import *
from dplot import dplot
from gmtpy import GMT, ScaleGuru, Ax, GridLayout, FrameLayout, inch
import pdb
from subprocess import *
import numpy as np
from numpy.random import random_integers
from minmax_mdl import minmax_mdl

def dissect_fname(fn):
    """
    get latitude, longitude and wave type (Love or Rayleigh) from
    filename.
    """
    lat = float(fn.split('_')[-2])
    lon = float(fn.split('_')[-1].split('.')[0] + '.' + fn.split('_')[-1].split('.')[1])
    # wtype = os.path.basename(fn).split('_')[1]
    # if wtype == 'love': wtype = 'L'
    # if wtype == 'rayleigh': wtype = 'R'
    return lat, lon

def make_widget():
    fig_width_pt = 483.7  # Get this from LaTeX using \showthe\columnwidth
    # fig_width_pt = 600.  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width / golden_mean  # height in inches
    fig_size = (fig_width, fig_height)
    fig_height = fig_height / 1.5
    fig_size = (fig_width, fig_height)

    masterlayout = GridLayout(2, 1)
    leftlayout = GridLayout(1, 3)
    rightlayout = GridLayout(1, 4)
    masterlayout.set_widget(0, 0, leftlayout)
    masterlayout.set_widget(1, 0, rightlayout)
    # # layout = GridLayout(2,4)
    widgets = []
    # ## dispersion curve widgets
    for _i in xrange(2):
        inner_layout = FrameLayout()
        leftlayout.set_widget(0, _i, inner_layout)
        widget = inner_layout.get_widget('center')
        widget.set_horizontal(0.5 * fig_size[0] * inch)
        widget.set_vertical(0.25 * fig_size[1] * inch)
        inner_layout.get_widget('top').set_vertical(0.05 * fig_size[1] * inch)
        inner_layout.get_widget('bottom').set_vertical(0.05 * fig_size[1] * inch)
        inner_layout.get_widget('right').set_horizontal(0.05 * fig_size[0] * inch)
        inner_layout.get_widget('left').set_horizontal(0.05 * fig_size[0] * inch)
        widgets.append(widget)

    # ## 1D profile widgets
    inner_layout = FrameLayout()
    leftlayout.set_widget(0, 2, inner_layout)
    widget = inner_layout.get_widget('center')
    widget.set_horizontal(0.5 * fig_size[0] * inch)
    widget.set_vertical(0.5 * fig_size[1] * inch)
    inner_layout.get_widget('top').set_vertical(0.1 * fig_size[1] * inch)
    inner_layout.get_widget('bottom').set_vertical(0.1 * fig_size[1] * inch)
    inner_layout.get_widget('left').set_horizontal(0.15 * fig_size[0] * inch)
    inner_layout.get_widget('right').set_horizontal(0.15 * fig_size[0] * inch)
    widgets.append(widget)

    # ## hist-widgets
    for _i in xrange(4):
        inner_layout = FrameLayout()
        rightlayout.set_widget(0, _i, inner_layout)
        widget = inner_layout.get_widget('center')
        widget.set_horizontal(0.5 * fig_size[0] * inch)
        widget.set_vertical(0.25 * fig_size[1] * inch)
        inner_layout.get_widget('top').set_vertical(0.05 * fig_size[1] * inch)
        inner_layout.get_widget('bottom').set_vertical(0.05 * fig_size[1] * inch)
        inner_layout.get_widget('right').set_horizontal(0.05 * fig_size[0] * inch)
        inner_layout.get_widget('left').set_horizontal(0.05 * fig_size[0] * inch)
        widgets.append(widget)
    return widgets, masterlayout

def get_disp(repfile, wtype='phase'):
    """
    read original dispersion curve from target file
    """
    fn = os.path.basename(repfile).replace('.report', '.target')
    # fn = fn.replace('cnipse_','')
    dirn = os.path.join(os.path.dirname(os.path.dirname(repfile)), 'targets')
    tfn = os.path.join(dirn, fn)
    p = []
    v = []
    err = []
    if wtype == 'phase':
        command = 'gptarget E -disp 0 %s' % tfn
        out = Popen(command, shell=True, stdout=PIPE).communicate()[0]
        d = out.split('\n')
    if wtype == "group":
        command = 'gptarget E -disp 1 %s' % tfn
        out = Popen(command, shell=True, stdout=PIPE).communicate()[0]
        d = out.split('\n')

    if wtype == "grouponly":
        command = 'gptarget E -disp 0 %s' % tfn
        out = Popen(command, shell=True, stdout=PIPE).communicate()[0]
        d = out.split('\n')

    if wtype == "joint":
        command = 'gptarget E -disp 2 %s' % tfn
        out = Popen(command, shell=True, stdout=PIPE).communicate()[0]
        d = out.split('\n')

    for _i in xrange(1, len(d)):
        try:
            a = d[_i].split()
            s = float(a[1]) * 1000.
            ds = float(a[2]) * 1000.
            p.append(1. / float(a[0]))
            v.append(1. / s)
            err.append(((1 / (s - ds) - 1 / s) + (1 / s - 1 / (s + ds)) / 2.))
        except Exception, e:
            print e
            continue
    return p, v, err

def get_best_model(repfile, dreso, dmax=40., dmin=0.):
    tmp = tempfile.mktemp()
    # ## read in models from dinver-output
    os.system('gpdcreport -best 1 %s >%s' % (repfile, tmp))
    mdl = loadtxt(tmp, skiprows=3)
    dh = 0
    nmdl = []
    nlayer = mdl.shape[0]
    for i in xrange(nlayer):
        th = mdl[i, 0] / 1000.
        vs = mdl[i, 2] / 1000.
        dold = dh
        nmdl.append([dold, vs])
        dh += th
        nmdl.append([dh, vs])
    a = array(nmdl)
    a[-1, 0] = dmax
    depth = arange(dmin, dmax + dreso, dreso)
    nvs = interp(depth, a[:, 0], a[:, 1])
    return vstack((nvs, -depth)).T


def plot_rep(repfile, paramfile, wtype, pixfile, show=True, misfit=0.1,
             minmax=True):
    """
    Plot 1D velocity profile.
    """
    # ## read in models from dinver-output
    os.system('gpdcreport ' + repfile + ' >tmp.mdl')
    models, data, nd, ne, nrow, nlayer = get_models('tmp.mdl')
    os.remove('tmp.mdl')
    if misfit == 'all':
        misfit = data.max()
    elif misfit <= data.min():
        misfit = median(data)
    # ## calculate density
    sreso = 0.05  # horizontal resolution
    dreso = .5  # vertical resolution
    grid1, grid2, x, y, smean, dmean = dplot(models, data, nd, ne, dreso=dreso,
                                             sreso=sreso, mf=misfit, nlayer=nlayer, dmax=120)
    # ## write output into temporary files that can be plotted by gmt
    xyz = tempfile.mktemp()
    xyz2 = tempfile.mktemp()
    grd = tempfile.mktemp()
    grdcpt = tempfile.mktemp()
    fout = pixfile
    f = open(xyz, 'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if grid1[ii, jj] > 0.0:
                print >> f, x[jj], -y[ii], grid1[ii, jj]
    f.close()
    f = open(xyz2, 'w')
    for ii in range(len(y)):
        for jj in range(len(x)):
            if grid2[ii, jj] > 0:
                print >> f, x[jj], -y[ii], '0.5'
    f.close()
    lat, lon = dissect_fname(repfile)
    if wtype == 'love':
        wt = 'L'
    if wtype == 'rayleigh':
        wt = 'R'
    # print repfile, lat, lon
    # ## gmt plot
    rng = '1/5/-120/0'
    step = '%f/%f' % (sreso, dreso)
    anot = int(grid1.max() / 1000.) * 1000 / 2.
    tick = anot / 2
    gmt = GMT(config={'ANOT_FONT_SIZE':8, 'LABEL_FONT_SIZE':10,
                      'ANNOT_OFFSET_SECONDARY':'0.1c',
                      'ANNOT_OFFSET_PRIMARY':'0.1c',
                      'LABEL_OFFSET':'0.1c',
                      'FRAME_PEN':'.5p'})
    widgets, layout = make_widget()
    if 1:
        widget = widgets[2]
        gmt.xyz2grd(xyz, G=grd, R=rng, I='%f/%f' % (sreso, dreso), out_discard=True)
        gmt.grd2cpt(grd, L='0/%d' % 3000, C="wysiwyg", D=True, Z=True, out_filename=grdcpt)
        gmt.psbasemap(R=True, B='a1f.5:S-velocity [km/s]:/a10f5:Depth [km]:WnSe', *widget.XYJ())
        gmt.psmask(xyz2, R=True, T=True, I='%f/%f' % (sreso, dreso), G='lightgray', *widget.XYJ())
        gmt.grdimage(grd, R=True, Q=True, C=grdcpt, *widget.XYJ())
        gmt.psxy(R=True, B=True, W='3,black', in_columns=[smean, -dmean], *widget.XYJ())
        bmdl = get_best_model(repfile, dreso, dmax=120., dmin=0.)
        gmt.psxy(R=True, B=True, W='3,black,-', in_rows=bmdl, *widget.XYJ())
        gmt.psscale(widget.XYJ()[0], widget.XYJ()[1], C=grdcpt, D='1.7c/1.5c/2.5c/.4ch', B='a%df%d10:No. of models:/::' % (1000, 500))
        txtstr = "1.2 -7. 8 0 1 LT lat = %3.2f" % lat
        gmt.pstext(R=rng, G='0/0/0', N=True, in_string=txtstr, *widget.XYJ())
        txtstr = "1.2 -12. 8 0 1 LT lon = %3.2f" % lon
        gmt.pstext(R=True, G='0/0/0', N=True, in_string=txtstr, *widget.XYJ())
        # plot the minimum and maximum model that is allowed
        # by the parameterisation
        if minmax:
            mdl_min, mdl_max = minmax_mdl(paramfile)
            gmt.psxy(R=True, B=True, W='3,red,-', in_rows=mdl_min, *widget.XYJ())
            gmt.psxy(R=True, B=True, W='3,red,-', in_rows=mdl_max, *widget.XYJ())
    if 1:
        p, v, e = get_disp(repfile, wtype='phase')
        plot_disp(gmt, repfile, widgets[1], p, v, e, 'phase', ptype='p', wtype=wt, mode=0, misfit=misfit)
        p, v, e = get_disp(repfile, wtype='group')
        plot_disp(gmt, repfile, widgets[0], p, v, e, 'group', ptype='g', wtype=wt, mode=0, misfit=misfit)
    if 0:
        p, v, e = get_disp(repfile, wtype='grouponly')
        plot_disp(gmt, repfile, widgets[0], p, v, e, 'group', ptype='g', wtype=wt, mode=0, misfit=misfit)
    if 0:
        plot_hist(gmt, widgets, models, nd, ne)
    gmt.save(fout)
    meanfout = fout.replace('.eps', '_mean.txt')
    savetxt(meanfout, vstack((-dmean, smean)).T)
    if show:
        os.system('gv %(fout)s&' % vars())


def plot_disp(gmt, mdl, widget, p, v, e, text, ptype='p', wtype='R', mode=0, misfit=0.1):
    """
    plot dispersion curves from inversion
    """
    # fn = tempfile.mktemp()
    fn = '/tmp/tmp.txt'
    if not os.path.isfile(mdl):
        print "%s does not exist"
        sys.exit(0)
    os.system('gpdcreport -best 1000 -%s%s %d %s >%s' % (ptype, wtype, mode, mdl, fn))
    nmdl = 100
    f = open(fn)
    nemax = int(f.readline().split()[1])
    curpos = f.tell()
    f.readline()
    nrow = 0
    while True:
        line = f.readline()
        if line.find('#') != -1: break
        a = line.split()
        nrow += 1
    f.seek(curpos)
    models = zeros((nemax, nrow, 2))
    data = zeros(nemax)
    ne = 0
    while True:
        line = f.readline()
        if not line: break
        if line.find('value') != -1:
            a = line.split()
            mf = float(a[9].split('=')[1])
            if mf < misfit:
                for _i in xrange(nrow):
                    line = f.readline()
                    a = line.split()
                    models[ne, _i, 0] = 1. / float(a[0])
                    models[ne, _i, 1] = 1. / (float(a[1]) * 1000.)
                data[ne] = mf
                ne += 1
    f.close()

    # try:
    mfmax = data[0:ne].max()
    mfmin = data[0:ne].min()
    # mfmin = 0.1
    # mfmax = 0.2
    step = (mfmax - mfmin) / 5.
    rng = '3/90/1.5/5'
    sca = 'X10.2c/4c'
    anot = 'a5f1:Period [s]:/a1f.5:Velocity [km/s]:WSne'
    cptf = tempfile.mktemp('tmp.cpt')
    x, y, j = widget.XYJ()
    x = float(x.split('-Xa')[1].split('p')[0])
    y = float(y.split('-Ya')[1].split('p')[0])
    j = float(j.split('-JX')[1].split('p')[0])
    print mfmin, mfmax, step
    scld = '%fp/%fp/%fp/.25ch' % (x + 50, y + 120, 0.3 * j)
    sclb = 'a%ff%f::/:misfit:' % (.001, .0005)
    zvals = tempfile.mktemp()
    savetxt(zvals, data)
    gmt.makecpt(C='wysiwyg', D=True, T='%f/%f/%f' % (mfmin, mfmax, step), I=True, out_filename=cptf)
    # gmt.makecpt(C='wysiwyg',D=True,T=zvals,I=True,out_filename=cptf)
    gmt.pstext(R=rng, G='0/0/0', in_string="%f %f 10 0 1 CM %s" % (80, 4.5, text), *widget.XYJ())
    for _i in reversed(sorted(data[random_integers(0, ne - 1, nmdl)].argsort())):
        dd = ['> ', '-Z %f' % data[_i]]
        dd = reshape(dd, (1, 2))
        gmt.psxy(R=True, C=cptf, M=True, B=anot, W='2', in_rows=append(dd, models[_i], axis=0), *widget.XYJ())
    # except:
    # pdb.set_trace()
    gmt.psxy(R=True, B=True, E='y/.5p', S='c.03', G='black', in_rows=vstack((vstack((p, v)), e)).T, *widget.XYJ())
    gmt.psscale(C=cptf, D=scld, B=sclb)


if __name__ == '__main__':
    savedir = '/home/behry/uni/data/results_pgc'
    result = os.path.join(savedir, 'reports/canada_48.000000_-105.000000.report')
    paramf = os.path.join(savedir, 'param/Lith50_Canada_48.000000_-105.000000.param')
    fout = os.path.join(savedir, 'plots', os.path.basename(result).replace('.report', '.eps'))
    misfit = 0.03
    wtype = 'rayleigh'
    plot_rep(result, paramf, wtype, fout, show=False, misfit='all')
    if 0:
        foutpdf = fout.replace('.ps', '.pdf')
        os.system('epstopdf --outfile=%s %s' % (foutpdf, fout))
        os.system("pdfcrop -margins '0. 0. 0. 0.' %s %s" % (foutpdf, foutpdf.replace('.pdf', '_crop.pdf')))

