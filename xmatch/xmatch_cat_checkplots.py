def xmatch_cat_checkplots(table1=None,
                          table2=None,
                          colnames_radec1=['ra', 'dec'],
                          colnames_radec2=['ra', 'dec'],
                          units_radec1=['degree', 'degree'],
                          units_radec2=['degree', 'degree'],
                          rmax=10.0,
                          rmax2=None,
                          plotfile_label=None,
                          plotfile_prefix=None,
                          suptitle=None,
                          title=None,
                          debug=False,
                          verbose=False):

    from xmatch import xmatch_checkplot1
    from xmatch import xmatch_checkplot2

    now = time.localtime(time.time())
    datestamp = time.strftime("%Y%m%d", now)
    function_name = inspect.stack()[0][3]

    lineno = str(inspect.stack()[0][2])
    #print(mk_timestamp(), function_name, lineno + ':')
    print(function_name + '.saveplot:', saveplot)
    print(function_name + '.plotfile:', plotfile)
    print(function_name + '.prefix:  ', plotfile_prefix)


    ra1 = table1[colnames_radec1[0]]
    dec1 = table1[colnames_radec1[1]]

    ra2 = table2[colnames_radec2[0]]
    dec2 = table2[colnames_radec2[1]]

    if plotfile_label is None:
        plotfile_label = ''

    if plotfile_prefix is None:
        plotfile_prefix = ''

    if plotfile_prefix is not None:
        plotfile_prefix = plotfile_prefix

    # suptitle = plotfile_label + 'nthN:' + str(nthneighbor)
    # suptitle = plotfile_label
    plotfile = (plotfile_prefix + 'xmatch_cat' +
                '_checkplot_1.png')

    # forked from Sophie Reed
    xmatch_checkplot1(
        ra1, dec1,
        ra2, dec2,
        width=rmax,
        gtype='square',
        saveplot=True,
        plotfile=plotfile,
        suptitle=suptitle)

    plotfile = (plotfile_prefix + 'xmatch_cat' +
                '_checkplot_2.png')

    # forked from Chris Desira
    xmatch_checkplot2(
                  ra1, dec1,
                  ra2, dec2,
                  width=rmax,
                  binsize=rmax/100,
                  gtype='square',
                  saveplot=True,
                  plotfile=plotfile,
                  suptitle=suptitle)

    return
