def xmatch_checkplots(ra1=None, dec1=None,
                      ra2=None, dec2=None,
                      dra=None, ddec=None, dr=None,
                      table1=None,
                      table2=None,
                      colnames_radec1=['ra', 'dec'],
                      colnames_radec2=['ra', 'dec'],
                      units_radec1=['degree', 'degree'],
                      units_radec2=['degree', 'degree'],
                      plotfile_label=None,
                      plotfile_prefix=None,
                      suptitle=None,
                      title=None,
                      rmax=10.0,
                      rmax2=None,
                      showplot=True,
                      saveplot=True,
                      datestamp=False,
                      verbose=False,
                      debug=False):
    """
    RA, Dec crossmatch validation plots based on code from Chris Desira and
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

    from xmatch import xmatch_checkplot1
    from xmatch import xmatch_checkplot2

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

    if plotfile_prefix is not None:
        plotfile_prefix = plotfile_prefix + '_'

    # suptitle = plotfile_label + 'nthN:' + str(nthneighbor)
    # suptitle = plotfile_label
    plotfile = (plotfile_prefix + 'xmatch_' + plotfile_label +
                '_checkplot_1.png')

    # forked from Sophie Reed
    xmatch_checkplot1(
        ra1, dec1,
        ra2, dec2,
        width=rmax,
        gtype='square',
        showplot=showplot,
        saveplot=saveplot,
        plotfile=plotfile,
        suptitle=suptitle)

    plotfile = (plotfile_prefix + 'xmatch_' + plotfile_label +
                '_checkplot_2.png')

    # forked from Chris Desira
    xmatch_checkplot2(
                  ra1, dec1,
                  ra2, dec2,
                  width=rmax,
                  binsize=rmax/100,
                  gtype='square',
                  saveplot=saveplot,
                  showplot=showplot,
                  plotfile=plotfile,
                  suptitle=suptitle)

    return
