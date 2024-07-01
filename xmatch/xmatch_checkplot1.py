from __future__ import (division, print_function)
#  Forked from Sophie Reed's version on 20160319

import os
import sys
import time

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std

from rgm_util import get_githash

def xmatch_checkplot1(
        dr=None, dra=None, ddec=None,
        ra1=None, dec1=None,
        ra2=None, dec2=None,
        units_radec1=['degree', 'degree'],
        units_radec2=['degree', 'degree'],
        figsize=(7.0, 7.0),
        hexbin=False,
        logfrequency=False,
        width=10.0,
        gtype="all",
        add_plotid=False,
        prefix=None,
        saveplot=True,
        showplots=True,
        title=None,
        suptitle=None,
        plotfile=None,
        plotfile_prefix=None,
        githash=False):

    """ Makes checkplot for catalogue xmatch results

    Forked from Sophie Reed's version on 20160319

    uses hist2d; a point based option would be useful

    Need to plot the circle that incribes the square

    Plot can either be square, the square inscribes the circle.
    Or all which has all the points in the matching circle.
    Square make the histograms more comparable.

    Compares RA_main and DEC_main columns with RA and Dec columns in the
    format output by the matching codes. Eg. RA_ + survey.

    Width needs to be in arcsecs

    """
    import math
    import time
    import inspect

    import numpy as np

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LogNorm

    # import stats
    # import plotid

    now = time.localtime(time.time())
    datestamp = time.strftime("%Y%m%d", now)
    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    function_name = inspect.stack()[0][3]

    lineno = str(inspect.stack()[0][2])
    print(timestamp, function_name, lineno + ':')
    print(function_name + '.saveplot:', saveplot)
    print(function_name + '.plotfile:', plotfile)
    print(function_name + '.prefix:  ', plotfile_prefix)

    if ra1 is not None and ra1.unit is None:
        ra1 = ra1 * u.deg
        print('len(ra1):', len(ra1), ra1.unit)
    if dec1 is not None and dec1.unit is None:
        dec1 = dec1 * u.deg
        print('len(dec1):', len(dec1), dec1.unit)

    if ra2 is not None and ra2.unit is None:
        ra2 = ra2 * u.deg
        print('len(ra2):', len(ra2), ra2.unit)
    if dec2 is not None and dec2.unit is None:
        dec2 = dec2 * u.deg
        print('len(dec2):', len(dec2), dec2.unit)

    if dra is not None:
        dra = dra * u.deg
        print('len(dra):',len(dra), dra.unit)
    if ddec is not None:
        ddec = ddec * u.deg
        print('len(ddec):', len(ddec), ddec.unit)
    if dr is not None:
        dr = dr * u.deg
        print('len(dr):',len(dr), dr.unit)

    print('Plot width:', width)
    print()

    if ra1 is not None and ra2 is not None:
        # compute Delta RA and Delta Dec in arcsecs
        # ra, dec assumed in have astropy units of degrees
        skycoord1 = SkyCoord(ra1, dec1, unit=units_radec1)
        skycoord2 = SkyCoord(ra2, dec2, unit=units_radec2)
        print(f'skycoord1[0]: skycoord1[0]', len(skycoord1))
        print(f'skycoord2[0]: skycoord2[0]', len(skycoord2))

        dra, ddec = skycoord1.spherical_offsets_to(skycoord2)
        dr = skycoord1.separation(skycoord2)

        print()
        print('dRA range:', np.min(dra), np.max(dra))
        print('dDec range:', np.min(ddec), np.max(ddec))
        print('dR range:', np.min(dr), np.max(dr))

        # convert to arcsecs
        print('convert to arcsec')
        dra = dra.arcsec
        ddec = ddec.arcsec
        dr = dr.arcsec

        print()
        print('dRA range:', np.min(dra), np.max(dra))
        print('dDec range:', np.min(ddec), np.max(ddec))
        print('dR range:', np.min(dr), np.max(dr))

    #width = width * u.arcsec
    ndata = len(dr)
    rmax = width
    print(len(dra))

    # convert offsets to arc seconds
    #print()
    #print('dra.unit:', dra.unit)
    #print('ddec.unit:', ddec.unit)
    #print('dr.unit:', dr.unit)
    #print('rmax.unit:', rmax.unit)
    #print('width.unit:', width.unit)

    # convert to arcsecs
    #try:
    #    dra = dra.arcsecond
    #except:
    #    dra = dra.to(u.arcsec)

    #try:
    #    ddec = ddec.arcsec
    #except:
    #    ddec = ddec.to(u.arcsec)

    #try:
    #    dr = dr.arcsec
    #except:
    #    dr = dr.to(u.arcsec)

    #try:
    #    rmax = rmax.arcsec
    #except:
    #    rmax = rmax.to(u.arcsec)

    #try:
    #    width = width.arcsec
    #except:
    #    width= width.to(u.arcsec)

    #print()
    #print('dra.unit:', dra.unit)
    #print('ddec.unit:', ddec.unit)
    #print('dr.unit:', dr.unit)
    #print('rmax.unit:', rmax.unit)
    #print('width.unit:', width.unit)


    print()
    print('dRA range:', np.min(dra), np.max(dra))
    print('dDec range:', np.min(ddec), np.max(ddec))
    print('dR range:', np.min(dr), np.max(dr))

    # radial cicular limits and square RA, Dec limits
    print('rmax:', rmax)
    print('width:', width)
    print()
    #print('rmax.unit:', rmax.unit)
    #print('dr.unit:', dr.unit)
    itest_dr = dr < rmax
    # itest_dr = np.where(dr < rmax)
    itest_dradec = ((dra > -1.0 * width) &
                    (dra < width) &
                    (ddec > -1.0 * width) &
                    (ddec < width))

    print('Number within match radius:', len(itest_dr),
          len(dr[itest_dr]), rmax)

    ndata_halfrmax = len(dr[np.where(dr < (rmax*0.5))])
    ndata_rmax = len(dr[np.where(dr < rmax)])
    ndata_2rmax = len(dr[np.where(dr < (2*rmax))])
    ndata_4rmax = len(dr[np.where(dr < (4*rmax))])

    print('Number of match with r < 0.5*rmax:', ndata_halfrmax)
    print('Number of match with r < rmax:', ndata_rmax,
          ndata_rmax - ndata_halfrmax)
    print('Number of match with r < 2.0*rmax:', ndata_2rmax,
                    ndata_2rmax - ndata_rmax)
    print('Number of match with r < 4.0*rmax:', ndata_4rmax,
          ndata_4rmax - ndata_2rmax)

    ndata_dradec_max = len(dr[np.where(
        (np.abs(dra) < rmax) & (np.abs(ddec) < rmax))])
    print('Number of match with radec < abs(rmax):', ndata_dradec_max)

    dra = dra[itest_dradec]
    ddec = ddec[itest_dradec]

    dr = dr[itest_dr]

    ra_med = np.median(dra)
    dec_med = np.median(ddec)
    ra_mad_std = mad_std(dra)
    dec_mad_std = mad_std(ddec)

    # error on mean/median
    ra_median_error = ra_mad_std / math.sqrt(len(dr))
    dec_median_error = dec_mad_std / math.sqrt(len(dr))

    print("Number of matchs", len(dra))
    print("RA median offset", ra_med)
    print("Dec median offset", dec_mad_std)
    print("RA Sigma(MAD)", ra_mad_std)
    print("Dec Sigma(MAD)", dec_mad_std)

    print("RA median error [sqrt(n)]", ra_median_error,
          "Dec median error [sqrt[n]", dec_median_error)

    print("dRA range:", np.min(dra), np.max(dra))
    print("dDec range:", np.min(ddec), np.max(ddec))

    xlimits = np.asarray([-1.0 * width, width])
    ylimits = np.asarray([-1.0 * width, width])
    limits = np.asarray([xlimits, ylimits])
    print(xlimits[0], xlimits[1])
    print(dra.dtype)
    print(dra.shape)
    print(xlimits.dtype)
    print(xlimits.shape)
    # itest = (xs > xlimits[0] & xs < xlimits[1])
    # xs = xs[itest]
    # itest = (ys > ylimits[0] & ys < ylimits[1])
    # ys = ys[itest]

    print('limits:', limits)
    print('width:', width)

    # GridSpec(nrows, ncols, figure=None,
    #          left=None, bottom=None, right=None, top=None,
    #          wspace=None, hspace=None,
    #          width_ratios=None,height_ratios=None)
    gs = gridspec.GridSpec(2, 2,
                           hspace=0.05, wspace=0.05,
                           width_ratios=[2, 1], height_ratios=[1, 2])

    fig = plt.figure(figsize=figsize)

    # Delta RA Histogram
    ax1 = plt.subplot(gs[0])
    ax1.hist(dra, bins=100,
             color="r", range=xlimits)
    plt.axvline(0.0, linestyle='dashed')
    ax1.set_xlim(xlimits)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_ylabel("Number")

    # Delta RA versus Delta Dec distribution
    ax2 = plt.subplot(gs[2])
    # ax2.plot(xs, ys, "k+")
    if len(dra) > 1000:
        xdata = dra
        ydata = ddec
        if not hexbin:
            cmap = 'cubehelix' # Dave Green
            # cmap = 'binary'
            plt.hist2d(xdata, ydata,
                       bins=100,
                       cmap=cmap,
                   norm=LogNorm(),
                   range=limits)
            plt.grid('true')
        if hexbin:
            plt.hexbin(xdata, ydata,
                       gridsize=100, bins='log',
                       cmap='cubehelix', mincnt=1.0)

    else:
        plt.plot(dra, ddec, "k.", ms=4)

    plt.axvline(0.0, linestyle='dashed')
    plt.axhline(0.0, linestyle='dashed')
    ax2.set_ylim(-1*width, width)
    ax2.set_xlim(-1*width, width)
    ax2.set_xlabel('Delta RA (arcsec)')
    ax2.set_ylabel('Delta Dec (arcsec')

    #labels1 = ax2.get_xticks()
    #ax2.set_xticklabels(labels1, rotation=270)

    # suptitle covers all all the subplots
    if suptitle is None:
        fig.suptitle("Number of sources in XMatch: " +
                     str(ndata), fontsize='small')

    if suptitle is not None:
        fig.suptitle(suptitle + ': ' + str(ndata), fontsize='medium')

    # Delta Dec
    ax3 = plt.subplot(gs[3])
    print('limits:', limits)
    ax3.hist(ddec, bins=100, orientation="horizontal", color="r",
        range=ylimits)

    ax3.set_ylim(ylimits)
    # ax3.set_xlabel("Number")
    plt.axhline(0.0, linestyle='dashed')
    ax3.axes.get_yaxis().set_visible(False)
    labels2 = ax3.get_xticks()
    ax3.set_xticklabels(labels2, rotation=270)

    ax4 = plt.subplot(gs[1])
    x0 = 0.0
    fontsize = 'small'
    fontsize = 'medium'
    ax4.annotate("Number of matchs: " +
                 str(len(dra)), xy=(x0, 0.0), size=fontsize)
    ax4.annotate("No. of sources(1): " +
                 str(len(ra1)),
                 xy=(x0, 0.20), size=fontsize)
    ax4.annotate("No. of sources(2): " +
                 str(len(ra2)),
                 xy=(x0, 0.10), size=fontsize)

    ax4.annotate("Median dRA: {0:.4f}".format(ra_med) +
                 '"', xy=(x0, 0.90), size=fontsize)
    ax4.annotate("Median dDEC: {0:.4f}".format(dec_med) +
                 '"', xy=(x0, 0.80), size=fontsize)
    ax4.annotate("dRA sigma MAD: {0:.4f}".format(ra_mad_std) +
                 '"', xy=(x0, 0.70), size=fontsize)
    ax4.annotate("dDEC sigma MAD: {0:.4f}".format(dec_mad_std) +
                 '"', xy=(x0, 0.60), size=fontsize)
    ax4.annotate("dRA median error: {0:.4f}".
                 format(ra_median_error) + '"',
                 xy=(x0, 0.50), size=fontsize)
    ax4.annotate("dDEC median error: {0:.4f}".
                 format(dec_median_error) + '"',
                 xy=(x0, 0.40), size=fontsize)

    ax4.axes.get_xaxis().set_visible(False)
    ax4.axes.get_yaxis().set_visible(False)
    ax4.axis('off')
    ax4.set_axis_off()


    if saveplot:
        lineno = str(inspect.stack()[0][2])
        #print(mk_timestamp(), function_name, lineno)
        print('plotfile:', plotfile)
        print('plotfile_prefix:', plotfile_prefix)
        if add_plotid:
            # make room for the plotid on right edge
            fig.subplots_adjust(right=0.95)
            # plotid()

        if plotfile is None:
            plotfile = 'match'
        if plotfile_prefix is not None and plotfile is None:
            plotfile = plotfile_prefix + '_xmatch_' + datestamp + '.png'
        if plotfile_prefix is None and plotfile is None:
            plotfile = 'xmatch_' + datestamp + '.png'

        # Add footnote; bbox is a dictionary
        footnote = timestamp + ': ' + \
        os.path.basename(__file__) + ': ' + plotfile
        plt.figtext(0.01, 0.01, footnote,
                    ha="left", fontsize=8,
                    bbox={"facecolor":'none', 'edgecolor':'none',
                          "alpha":0.5, "pad":2})


        if githash:
            git_hash = get_githash()
            text = git_hash
            xtext = 0.01
            ytext = 0.97
            transform = plt.gcf().transFigure
            color = 'b'
            fontsize = 8
            plt.figtext(xtext, ytext,
                        text,
                        transform=transform,
                        rotation=90.0,
                        color=color,
                        backgroundcolor='w',
                        fontsize=fontsize,
                        weight='ultralight',
                        horizontalalignment='left',
                        verticalalignment='bottom')



        print('Saving: ', plotfile)
        plt.savefig(plotfile)


    if showplots:
        if not saveplot:
            footnote = timestamp + ': ' + \
                os.path.basename(__file__)
            plt.figtext(0.01, 0.01, footnote,
                        ha="left", fontsize=8,
                        bbox={"facecolor":'none',
                              "edgecolor":'none',
                              "alpha":0.5, "pad":2})
        plt.show()


    return ra_med, dec_med
