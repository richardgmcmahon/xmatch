def xmatch_cat_checkplots(table1=None,
                          table2=None,
                          colnames_radec1=['ra', 'dec'],
                          colnames_radec2=['ra', 'dec'],
                          units_radec1=['degree', 'degree'],
                          units_radec2=['degree', 'degree'],
                          plotfile_label=None,
                          suptitle=None,
                          title=None,
                          rmax=10.0,
                          rmax2=None,
                          debug=False,
                          verbose=False):

    from xmatch import xmatch_checkplot1
    from xmatch import xmatch_checkplot2

    ra1 = table1[colnames_radec1[0]]
    dec1 = table1[colnames_radec1[1]]

    ra2 = table2[colnames_radec2[0]]
    dec2 = table2[colnames_radec2[1]]

    if plotfile_label is None:
        plotfile_label = ''

    # suptitle = plotfile_label + 'nthN:' + str(nthneighbor)
    # suptitle = plotfile_label
    plotfile = 'xmatch_cat' + plotfile_label + '_checkplot_1.png'

    # forked from Sophie Reed
    xmatch_checkplot1(
        ra1, dec1,
        ra2, dec2,
        width=rmax,
        gtype='square',
        saveplot=True,
        plotfile=plotfile,
        suptitle=suptitle)

    plotfile = 'xmatch_cat' + plotfile_label + '_checkplot_2.png'

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
