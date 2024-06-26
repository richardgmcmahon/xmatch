def xmatch_checkplots(dr=None, dra=None, ddec=None,
                      ra1=None, dec1=None,
                      ra2=None, dec2=None,
                      table1=None,
                      table2=None,
                      colnames_radec1=['ra', 'dec'],
                      colnames_radec2=['ra', 'dec'],
                      units_radec1=['degree', 'degree'],
                      units_radec2=['degree', 'degree'],
                      rmax=10.0,
                      rmax2=None,
                      plotfile_label=None,
                      plotfile_prefix=None,
                      plotfile_suffix=None,
                      suptitle=None,
                      title=None,
                      showplots=True,
                      saveplot=True,
                      datestamp=False,
                      hexbin=False,
                      githash=False,
                      debug=False,
                      verbose=False):
    """ RA, Dec crossmatch validation plots based on code from Chris Desira and
    Sophie Reed

    Docstring follows the Pandas convention
    https://pandas.pydata.org/docs/development/contributing_docstring.html
    which is based on: https://numpydoc.readthedocs.io/en/latest/format.html


Parameters
    ----------
    ra1: real
        Right Ascension or Longitude in degrees for catalogue or table #1
    dec1: real
        Declination or Latitude in degrees for catalogue or table #1

    ra2: real
        Right Ascension or Longitude in degrees for catalogue or table #1
    dec2: real
        Declination or Latitude in degrees for catalogue or table #1

    Astropy units are supported so radians can be passed transparently

    Returns
    -------
    int
        The sum of ``num1`` and ``num2``.

    See Also
    --------
    subtract : Subtract one integer from another.

    Examples
    --------
    >>> add(2, 2)
    4
    >>> add(25, 0)
    25
    >>> add(10, -10)
    0


    """

    import os
    import sys
    import time
    import inspect

    from xmatch import xmatch_checkplot1
    from xmatch import xmatch_checkplot2

    from rgm_util import get_githash

    now = time.localtime(time.time())
    datestamp = time.strftime("%Y%m%d", now)
    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    function_name = inspect.stack()[0][3]

    lineno = str(inspect.stack()[0][2])
    print(timestamp, function_name, lineno + ':')
    print(function_name + '.saveplot:', saveplot)
    print(function_name + '.prefix:  ', plotfile_prefix)
    print(function_name + '.suffix:  ', plotfile_suffix)

    if dr is None:
        if ra1 is None:
            ra1 = table1[colnames_radec1[0]]
        if dec1 is None:
            dec1 = table1[colnames_radec1[1]]

        if ra2 is None:
            ra2 = table2[colnames_radec2[0]]
        if dec2 is None:
            dec2 = table2[colnames_radec2[1]]

    if plotfile_label is None:
        plotfile_label = ''

    if plotfile_prefix is None:
        plotfile_prefix = ''

    # suptitle = plotfile_label + 'nthN:' + str(nthneighbor)
    # suptitle = plotfile_label
    plotfile = plotfile_prefix + 'xmatch_checkplot_1.png'

    if plotfile_suffix is not None:
        plotfile = plotfile_prefix + 'xmatch_checkplot_1' + \
         plotfile_suffix + '.png'


    # forked from Sophie Reed
    title = plotfile_prefix
    suptitle = plotfile_prefix
    xmatch_checkplot1(
        dr=dr, dra=dra, ddec=ddec,
        ra1=ra1, dec1=dec1,
        ra2=ra2, dec2=dec2,
        width=rmax,
        hexbin=hexbin,
        gtype='square',
        saveplot=True,
        plotfile=plotfile,
        title=title,
        suptitle=suptitle,
        showplots=showplots,
        githash=githash)


    plotfile = (plotfile_prefix + 'xmatch_checkplot_2.png')
    if plotfile_suffix is not None:
        plotfile = plotfile_prefix + 'xmatch_checkplot_2' + \
         plotfile_suffix + '.png'


    # forked from Chris Desira
    xmatch_checkplot2(
        dr=dr, dra=dra, ddec=ddec,
        ra1=ra1, dec1=dec1,
        ra2=ra2, dec2=dec2,
        hexbin=hexbin,
        width=rmax,
        binsize=rmax/100,
        gtype='square',
        saveplot=True,
        plotfile=plotfile,
        suptitle=suptitle,
        showplots=showplots,
        githash=githash)

    return
