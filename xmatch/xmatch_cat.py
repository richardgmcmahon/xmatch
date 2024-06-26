from __future__ import (division, print_function)

def xmatch_cat(ra1=None, dec1=None,
               ra2=None, dec2=None,
               table1=None, table2=None,
               radec1=None, radec2=None,
               colnames_radec1=['ra', 'dec'],
               colnames_radec2=['ra', 'dec'],
               units_radec1=['degree', 'degree'],
               units_radec2=['degree', 'degree'],
               nthneighbor=None,
               multimatch=False,
               seplimit=10.0,
               selfmatch=False,
               stats=True,
               debug=False,
               verbose=False,
               method=False):
    """RA, Dec nearest xmatch for two lists; returns pointers

    returns

    if not multimatch:
        return idx2, dr, dra, ddec

    if multimatch:
        return (idx1, idx2), dr, dra, ddec


    nearest match

    I am not sure why I adopted the zipped option; maybe it is for speed
    when using astropy. radec1, radec2 are not currently supported!

    input can be an astropy table or zipped radec as a list

    e.g.

    c = zip([1],[1])
    radec1 = zip(ra1 , dec1)


    radec1 = np.column_stack(ra1, dec1))


    Self match notes:


    Known issues:
    Some table will have wrong units e.g. None or radians instead of radian

    To change try something like:
    table[colname_RA].unit = 'radian'
    table[colname_Dec].unit = 'radian'


    """

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.table import Table, hstack
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import search_around_sky, match_coordinates_sky
    from astropy import units as u

    from astropy.stats import mad_std, median_absolute_deviation

    if verbose or debug:
        print('__file__:', __file__)
        print('__name__:', __name__)

    try:
        if 'filename' in table1.meta:
            print('table1.filename:', table1.meta['filename'])
    except:
        print("table1 has no metadata or table1.meta['filename']")

    if verbose or debug:
        print('colnames_radec1:', colnames_radec1)
        table1.info()

    # selfmatch does not need a 2nd table
    if not selfmatch:
        try:
            if 'filename' in table2.meta:
                print('table2.filename:', table2.meta['filename'])
        except:
            print("table2 has no metadata or table2.meta['filename']")

        if verbose or debug:
            print('colnames_radec2:', colnames_radec2)
            table2.info()

    if selfmatch:
        table2 = table1
        colnames_radec2 = colnames_radec1
        if nthneighbor is None:
            nthneighbor = 2

    if nthneighbor is None:
        nthneighbor = 1

    if ra1 is None:
        ra1 = table1[colnames_radec1[0]]
    if dec1 is None:
        dec1 = table1[colnames_radec1[1]]
    if verbose or debug:
        print('table1: ra; column name, units:',
              colnames_radec1[0], table1[colnames_radec1[0]].unit)
        print('table1: dec: column name, units:',
              colnames_radec1[1], table1[colnames_radec1[1]].unit)

    if ra2 is None:
        ra2 = table2[colnames_radec2[0]]
    if dec2 is None:
        dec2 = table2[colnames_radec2[1]]
    if verbose or debug:
        print('table2: ', colnames_radec2[0], table2[colnames_radec2[0]].unit)
        print('table2: ', colnames_radec2[1], table2[colnames_radec2[1]].unit)

    if stats or verbose or debug:
        print('RA1 range:', np.min(ra1), np.max(ra1))
        print('Dec1 range:', np.min(dec1), np.max(dec1))

        print('RA2 range:', np.min(ra2), np.max(ra2))
        print('Dec2 range:', np.min(dec2), np.max(dec2))

    print('units_radec1:', units_radec1)
    print('units_radec2:', units_radec2)

    if ra1.unit is None:
        ra1.unit = units_radec1[0]
    if dec1.unit is None:
        dec1.unit = units_radec1[1]

    if ra2.unit is None:
        ra2.unit = units_radec2[0]
    if dec1.unit is None:
        dec2.unit = units_radec2[1]

    print('ra1.unit:', ra1.unit)
    print('dec1.unit:', dec1.unit)
    print('ra2.unit:', ra2.unit)
    print('dec2.unit:', dec2.unit)

    if verbose or debug:
        print('Convert to SkyCoord', len(ra1), len(ra2))

    skycoord1 = SkyCoord(ra1, dec1, unit=units_radec1, frame='icrs')
    skycoord2 = SkyCoord(ra2, dec2, unit=units_radec2, frame='icrs')

    # idx is an integer array into the second cordinate array to get the
    # matched points for the second coordindate array.
    # Shape of idx matches the first coordinate array
    idx1 = []
    idx2 = []
    if verbose or debug:
        print('method:', method)
        print('multimatch:', multimatch)

    if not method:
        if not multimatch:
            idx2, d2d, d3d = \
                match_coordinates_sky(skycoord1,
                                      skycoord2,
                                      nthneighbor=nthneighbor)
        if multimatch:
            idx1, idx2, d2d, d3d = \
                search_around_sky(skycoord1,
                                  skycoord2,
                                  seplimit * u.arcsec)

    # alternative 'method' form
    if method:
        if not multimatch:
            idx2, d2d, d3d = \
                skycoord1.match_to_catalog_sky(skycoord2,
                                              nthneighbor=nthneighbor)

        if multimatch:
            idx1, idx2, d2d, d3d = \
                skycoord1.search_around_sky(skycoord2,
                                            seplimit * u.arcsec)

    # apply seplimit nearest match case
    #if not multimatch:
    #    itest = d2d < seplimit * u.arcsec
    #    print('Number of matchs within seplimit:',
    #          len(d2d[itest]), len(d2d), seplimit)
    #    idx2 = idx2[itest]
    #    d2d = d2d[itest]
    #    d3d = d3d[itest]

    # compute the separations and
    if verbose or debug:
        print('Compute separations')
    if not multimatch:
        separation = skycoord1.separation(skycoord2[idx2])
        dra, ddec = \
            skycoord1.spherical_offsets_to(skycoord2[idx2])


    if multimatch:
        separation = skycoord1[idx1].separation(skycoord2[idx2])
        dra, ddec = \
            skycoord1[idx1].spherical_offsets_to(skycoord2[idx2])


    if stats or verbose or debug:
        print('multimatch:', multimatch)
        print('seplimit:', seplimit)
        print('len(table1):', len(table1))
        print('len(table2):', len(table2))
        print()

        print('len(idx1):', len(idx1))
        print('len(idx2):', len(idx2))
        print('idxmatch range:', np.min(idx2), np.max(idx2))
        print()

        print('d2d min, max:', np.min(d2d), np.max(d2d))
        print('d2d min, max (arcsec):', np.min(d2d).arcsec, np.max(d2d).arcsec)
        print('d2d median:', np.median(d2d).arcsec)

        median_separation = np.median(separation).arcsec
        mad_std_separation = mad_std(separation.arcsec)
        print('dR range (arcsec):',
              np.min(separation.arcsec), np.max(separation.arcsec))
        print('dR mean, std (arcsec):',
              np.mean(separation).arcsec, np.std(separation).arcsec)
        print('dR  median, mad_std (arcsec):',
              median_separation, mad_std_separation)
        print()

        median_dra = np.median(dra).arcsec
        mad_std_dra = mad_std(dra.arcsec)
        print('dRA min, max:',
              np.min(dra).arcsec, np.max(dra).arcsec)
        print('dRA mean, std:',
              np.mean(dra).arcsec, np.std(dra).arcsec)
        print('dRA median, mad_std:',
              median_dra, mad_std_dra)
        print()

        median_ddec = np.median(ddec).arcsec
        mad_std_ddec = mad_std(ddec.arcsec)
        print('dDec min, max:',
              np.min(ddec).arcsec, np.max(ddec).arcsec)
        print('dDec mean, std:',
              np.mean(ddec).arcsec, np.std(ddec).arcsec)
        print('dDec median, mad_std:',
              median_ddec, mad_std_ddec)
        print()

    # convert to arcsec for convenience
    separation = separation.arcsec
    dr = d2d.arcsec
    dra = dra.arcsec
    ddec = ddec.arcsec

    # return dra, ddec, dr in arcsec
    # as a list or could be dict; check if scales from 10^3 -> 10^6 -> 10^9
    drplus = [dra, ddec, dr]

    if debug or verbose:
        print(len(idx2), len(dr))
        print(len(drplus), len(drplus[0]), len(drplus[1]), len(drplus[2]))

    # could add option to return dr, dra, ddec
    if not multimatch:
        return idx2, dr, dra, ddec

    if multimatch:
        return (idx1, idx2), dr, dra, ddec
