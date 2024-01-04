from __future__ import (division, print_function)
#  Forked from Sophie Reed's version on 20160319
import time

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std


def xmatch_checkplot1(ra1, dec1,
                      ra2, dec2,
                      units_radec1=['degree', 'degree'],
                      units_radec2=['degree', 'degree'],
                      figsize = (7.0, 7.0),
                      width=10.0,
                      gtype="all",
                      add_plotid=True,
                      prefix=None,
                      saveplot=True,
                      showplot=True,
                      plotfile=None,
                      plotfile_prefix=None,
                      title=None,
                      suptitle=None):

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
    function_name = inspect.stack()[0][3]

    lineno = str(inspect.stack()[0][2])
    #print(mk_timestamp(), function_name, lineno + ':')
    print(function_name + '.saveplot:', saveplot)
    print(function_name + '.plotfile:', plotfile)
    print(function_name + '.prefix:  ', plotfile_prefix)
    print(len(ra1), len(ra2))

    print('Plot width:', width)
    ndata = len(ra1)
    rmax = width

    # compute Delta RA and Delta Dec in arcsecs
    # ra, dec assumed in have astropy units of degrees
    skycoord1 = SkyCoord(ra1, dec1, unit=units_radec1)
    skycoord2 = SkyCoord(ra2, dec2, unit=units_radec2)
    print(skycoord1[0])
    print(skycoord2[0])

    dra, ddec = skycoord1.spherical_offsets_to(skycoord2)
    dr = skycoord1.separation(skycoord2)

    print(len(dra))
    print(dra[0])
    print(ddec[0])

    # convert offsets to arc seconds
    dra = dra.arcsecond
    ddec = ddec.arcsecond
    dr = dr.arcsecond

    print('dRA range:', np.min(dra), np.max(dra))
    print('dDec range:', np.min(ddec), np.max(ddec))
    print('dR range:', np.min(dr), np.max(dr))

    # radial cicular limits and square RA, Dec limits
    itest_dr = np.where(dr < rmax)
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
    ax1.hist(dra, bins=100, color="r", range=xlimits)
    plt.axvline(0.0, linestyle='dashed')
    ax1.set_xlim(xlimits)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_ylabel("Number")


    # Delta RA versus Delta Dec distribution
    ax2 = plt.subplot(gs[2])
    # ax2.plot(xs, ys, "k+")
    if len(dra) > 1000:
        plt.hist2d(dra, ddec, bins=100,
                   cmap="binary",
                   norm=LogNorm(),
                   range=limits)
        plt.grid('true')
    else:
        plt.plot(dra, ddec, "k.", ms=4)


    plt.axvline(0.0, linestyle='dashed')
    plt.axhline(0.0, linestyle='dashed')
    ax2.set_ylim(-1*width, width)
    ax2.set_xlim(-1*width, width)
    ax2.set_xlabel('Delta RA "')
    ax2.set_ylabel('Delta Dec "')


    #labels1 = ax2.get_xticks()
    #ax2.set_xticklabels(labels1, rotation=270)

    # Delta Dec
    if suptitle is None:
        fig.suptitle("Number of sources in XMatch: " +
                     str(ndata), fontsize='small')

    if suptitle is not None:
        fig.suptitle(suptitle + ': ' + str(ndata), fontsize='small')

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
            plotfile = plotfile_prefix + '_match_' + datestamp + '.png'
        if plotfile_prefix is None and plotfile is None:
            plotfile = 'match_' + datestamp + '.png'

        print('Saving: ', plotfile)
        plt.savefig(plotfile)

    if showplot:
        plt.show()

    return ra_med, dec_med
