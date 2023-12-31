def xmatch_checkplot2(ra1, dec1,
                      ra2, dec2,
                      units_radec1=['degree', 'degree'],
                      units_radec2=['degree', 'degree'],
                      width=10.0,
                      binsize=0.1,
                      markersize=1.0,
                      showplot=True,
                      saveplot=True,
                      plotfile='',
                      plotfile_prefix=None,
                      suptitle=None,
                      **kwargs):

    """
    only use suptitle since title near start
    before subplots causes alignment problems

    """

    import time
    import inspect

    import numpy as np

    import matplotlib.pyplot as plt

    from astropy import stats
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    #try:
    #    from mk_timestamp import mk_timestamp
    #from plotid import plotid


    now = time.localtime(time.time())
    datestamp = time.strftime("%Y%m%d", now)
    function_name = inspect.stack()[0][3]

    lineno = str(inspect.stack()[0][2])
    #print(mk_timestamp(), function_name, lineno + ':')
    print(function_name + '.saveplot:', saveplot)
    print(function_name + '.plotfile:', plotfile)
    print(function_name + '.prefix:  ', plotfile_prefix)
    print(len(ra1), len(ra2))

    rmax = width

    if suptitle is None:
        suptitle=''

    ndata_all = len(ra1)

    print('RA1 range:', np.min(ra1), np.max(ra1))
    print('Dec1 range:', np.min(dec1), np.max(dec1))
    print('RA2 range:', np.min(ra2), np.max(ra2))
    print('Dec2 range:', np.min(dec2), np.max(dec2))

    # compute Delta RA and Delta Dec in arcsecs
    # ra, dec assumed in have astropy units of degrees
    skycoord1 = SkyCoord(ra1, dec1, unit=units_radec1)
    skycoord2 = SkyCoord(ra2, dec2, unit=units_radec2)

    dra, ddec = skycoord1.spherical_offsets_to(skycoord2)
    dr = skycoord1.separation(skycoord2)
    print(len(dra))
    print(dra[0])
    print(ddec[0])

    dra = dra.arcsecond
    ddec = ddec.arcsecond
    dr = dr.arcsecond

    rmax = width

    itest = (np.abs(dra) < rmax) & (np.abs(ddec) < rmax)

    itest_dr = np.where(dr < rmax)
    itest_dradec = ((dra > -1.0 * width) &
                    (dra < width) &
                    (ddec > -1.0 * width) &
                    (ddec < width))

    dra = dra[itest_dradec]
    ddec = ddec[itest_dradec]
    dr = dr[itest_dr]

    dr_median = np.median(dr)
    dr_ndata = len(dr)
    dr_mad = stats.median_absolute_deviation(dr)
    dr_mad_std = stats.mad_std(dr)

    dra_median = np.median(dra)
    dra_ndata = len(dra)
    dra_mad = stats.median_absolute_deviation(dra)
    dra_mad_std = stats.mad_std(dra)

    ddec_median = np.median(ddec)
    ddec_ndata = len(ddec)
    ddec_mad = stats.median_absolute_deviation(ddec)
    ddec_mad_std = stats.mad_std(ddec)

    # error on mean/median
    dra_median_error = dra_mad_std / np.sqrt(len(dr))
    ddec_median_error = ddec_mad_std / np.sqrt(len(dr))


    fig = plt.figure(1, figsize=(10, 5))
    plt.suptitle(suptitle + ': '+ str(ndata_all))
    ax1=fig.add_subplot(1,2,1)

    xdata = dr
    ndata = len(xdata)
    bins = int(rmax/binsize)
    n, b, patches = ax1.hist(xdata, bins=bins,
                             color='green', alpha=0.5)

    bin_min = np.where(n == n.min())

    ax1.locator_params(axis='x', nbins=4)

    s04 = '# = %i'% ndata
    ax1.annotate(s04,(0.28,0.90) , xycoords = 'axes fraction',size=8)

    s01 = 'Median separation = %.4f' % dr_median
    ax1.annotate(s01,(0.28,0.85) , xycoords = 'axes fraction',size=8)

    ax1.set_xlabel('Pairwise separation (arcseconds)')
    ax1.set_ylabel('Frequency per bin')

    ax2 = fig.add_subplot(1,2,2, aspect='equal')

    alpha = 1.0
    ndata = len(dra)
    ax2.plot(dra,ddec,'oc',
             markersize=markersize,
             markeredgewidth=0.0,
             alpha=alpha) #0.5 smallest size

    ax2.axis([-1.0*rmax, rmax,-1.0*rmax, rmax])
    ax2.locator_params(axis='x',nbins=4)
    ax2.set_xlabel('Delta RA')
    ax2.set_ylabel('Delta Dec')

    s1 = 'xmatchs: %i' % ndata
    ax2.annotate(s1,(0.45,0.95) , xycoords = 'axes fraction',size=8)

    s2 = 'dRA Median = %.4f' % dra_median
    ax2.annotate(s2,(0.05,0.90) , xycoords = 'axes fraction',size=8)
    s3 = 'dDec Median = %.4f' % ddec_median
    ax2.annotate(s3,(0.55,0.90) , xycoords = 'axes fraction',size=8)

    s2 = 'dRA sigma_MAD = %.4f' % dra_mad_std
    ax2.annotate(s2,(0.05,0.85) , xycoords = 'axes fraction',size=8)
    s3 = 'dDec sigma_MAD = %.4f' % ddec_mad_std
    ax2.annotate(s3,(0.55,0.85) , xycoords = 'axes fraction',size=8)


    s2 = 'dRA Median err = %.4f' % dra_median_error
    ax2.annotate(s2,(0.05,0.80) , xycoords = 'axes fraction',size=8)
    s3 = 'dDec Median err = %.4f' % ddec_median_error
    ax2.annotate(s3,(0.55,0.80) , xycoords = 'axes fraction',size=8)


    fig.tight_layout()
    ax2.grid()

    fig.subplots_adjust(top=0.88)


    # make room for the plotid on right edge
    fig.subplots_adjust(right=0.95)
    #plotid()

    if plotfile != None:
        print('Saving plotfile:', plotfile)
        plt.savefig(plotfile)

    if ('save' in kwargs):
        path_to_save = str(kwargs['save'])
        plt.savefig(path_to_save, dpi=150)
    else:
        if showplot:
            plt.show()

    return
