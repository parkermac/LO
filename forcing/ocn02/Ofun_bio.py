"""
Functions to add biogeochemical fields to a clm file.
"""

import numpy as np
import matplotlib.path as mpth
import pandas as pd
import gsw

verbose = False

# These values are mostly generated by LPM/ic/ic_test2.py, with DIC and TA
#
# Values for sog filled in by hand.
#
# The info is roganized as a dict of dicts, with the outer keys being basins,
# and the inner keys being variables names (conforming to LO/obs standards).
# The values are tuples of mean values from bottle data over some number of
# years and for selected months, and the means are for z (above,below)
# some dividing z.
#
# These numbers are for years 2013-2017, and months 9-11, with dividing z = -25 m
basin_clim_dict = {
'sog': {
    'salt':(28,30),
    'temp':(12,9),
    'oxygen':(226,141),
    'NO3':(17,27),
    'TIC':(2029,2098), # copied from ps
    'alkalinity':(2106,2149), # copied from ps
    },
'jdf': {
    'salt':(31,32),
    'temp':(10,9),
    'oxygen':(193,136),
    'NO3':(20,27),
    'TIC':(2116,2213),
    'alkalinity':(2185,2247),
    },
'ps': {
    'salt':(29,30),
    'temp':(13,12),
    'oxygen':(223,177),
    'NO3':(16,21),
    'TIC':(2029,2098),
    'alkalinity':(2106,2149),
    },
'hc': {
    'salt':(28,30),
    'temp':(12,11),
    'oxygen':(192,129),
    'NO3':(17,25),
    'TIC':(2041,2124),
    'alkalinity':(2083,2143),
    },
}


def fill_polygons(fld, vn, G, z_rho, Ldir):
    """
    This takes a 4-D array [t,z,y,x] that is on its way to being an ocean_clm.nc field
    and overwrites values in the Salish Sea, by basin, to appropriate climatological
    values.
    
    The development code for this function is in LPM/ic.
    """
    # grid
    x = G['lon_rho']
    y = G['lat_rho']
    xy = np.concatenate((x.reshape(-1,1),y.reshape(-1,1)), axis=1)

    # polygons
    basin_list = ['sog', 'jdf', 'ps','hc']
    bb = 0
    for basin in basin_list:
        # polygon
        fnp = Ldir['LO'] / 'forcing' / Ldir['frc'] / 'polygons' / ('poly_'+basin+'.p')
        p = pd.read_pickle(fnp)
        xx = p.x.to_numpy()
        yy = p.y.to_numpy()
        xxyy = np.concatenate((xx.reshape(-1,1),yy.reshape(-1,1)), axis=1)
        path = mpth.Path(xxyy)
        
        # make Boolean array [y,x] of which are inside
        isin = path.contains_points(xy)
        isina = isin.reshape(x.shape)
        
        zdiv = -25
        vtop = basin_clim_dict[basin][vn][0]
        vbot = basin_clim_dict[basin][vn][1]
        zmask = z_rho >= zdiv
                
        if bb == 0:
            fld_temp = fld[0,:,:,:].copy()
            fld_temp[zmask] = vtop
            fld_temp[~zmask] = vbot
            fld_new = fld[0,:,:,:].copy()
        
        fld_new[:,isina] = fld_temp[:,isina]
        
        bb += 1
        
    nt = fld.shape[0]
    fld_new_full = fld_new[None,:] * np.ones((nt,1,1,1))
        
    return fld_new_full
    
        

def create_bio_var(salt, vn):
    # check salt limits
    if np.nanmax(salt) >= 35:
        print('Warning from Ofun_bio.create_bio_var(): salt > 35')
        print('Max salt value = %0.2f' % (np.nanmax(salt)))
        print('Action taken: limited salt in calculation to 34.99')
        salt[salt>=35] = 34.99
    
    if verbose:
        print('  -- adding ' + vn)
    if vn == 'NO3':
        # Salinity vs. NO3 [uM], Ryan McCabe 8/2015
        # NO3 = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.898)
        mm[ind] = 0
        bb[ind] = 0
        ind = ((salt >= 31.898) & (salt < 33.791))
        mm[ind] = 16.3958
        bb[ind] = -522.989
        ind = ((salt >= 33.791) & (salt < 34.202))
        mm[ind] = 29.6973
        bb[ind] = -972.4545
        ind = ((salt >= 34.202) & (salt < 34.482))
        mm[ind] = 8.0773
        bb[ind] = -233.0007
        ind = ((salt >= 34.482) & (salt < 35))
        mm[ind] = -28.6251
        bb[ind] = 1032.5686
        NO3 = mm*salt + bb
        # Set maximum NO3 to 45 microMolar (found at ~800m depth), based on
        # evidence from historical NO3 data in NODC World Ocean Database.
        NO3[NO3 > 45] = 45
        # Ensure that there are no negative values.
        NO3[NO3 < 0] = 0
        return NO3
    elif vn == 'NH4':
        NH4 = 0 * salt
        return NH4
    elif vn == 'oxygen':
        if np.nanmax(salt) > 36:
            print('Salt out of range for oxgen regression')
        # Salinity vs. oxygen [uM], Ryan McCabe 8/2015
        # oxygen = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 32.167)
        mm[ind] = 0
        bb[ind] = 300
        ind = ((salt >= 32.167) & (salt < 33.849))
        mm[ind] = -113.9481
        bb[ind] = 3965.3897
        ind = ((salt >= 33.849) & (salt < 34.131))
        mm[ind] = -278.3006
        bb[ind] = 9528.5742
        ind = ((salt >= 34.131) & (salt < 34.29))
        mm[ind] = -127.2707
        bb[ind] = 4373.7895
        ind = ((salt >= 34.29) & (salt < 34.478))
        mm[ind] = 34.7556
        bb[ind] = -1182.0779
        ind = ((salt >= 34.478) & (salt < 35))
        mm[ind] = 401.7916
        bb[ind] = -13836.8132
        oxygen = mm*salt + bb
        # Limit values.
        oxygen[oxygen > 450] = 450
        oxygen[oxygen < 0] = 0
        return oxygen
    elif vn == 'TIC':
        if np.nanmax(salt) > 36:
            print('Salt out of range for TIC regression')
        # Salinity vs. TIC [uM]
        # TIC = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.887)
        mm[ind] = 27.7967
        bb[ind] = 1112.2027
        ind = ((salt >= 31.887) & (salt < 33.926))
        mm[ind] = 147.002
        bb[ind] = -2688.8534
        ind = ((salt >= 33.926) & (salt < 34.197))
        mm[ind] = 352.9123
        bb[ind] = -9674.5448
        ind = ((salt >= 34.197) & (salt < 34.504))
        mm[ind] = 195.638
        bb[ind] = -4296.2223
        ind = ((salt >= 34.504) & (salt < 35))
        mm[ind] = -12.7457
        bb[ind] = 2893.77
        TIC = mm*salt + bb
        return TIC
    elif vn == 'alkalinity':
        mm = 0*salt
        bb = 0*salt
        if np.nanmax(salt) > 36:
            print('Salt out of range for alkalinity regression')
        # Salinity vs. alkalinity [uM]
        # alkalinity = mm*salt + bb
        ind = (salt < 31.477)
        mm[ind] = 37.0543
        bb[ind] = 1031.0726
        ind = ((salt >= 31.477) & (salt < 33.915))
        mm[ind] = 48.5821
        bb[ind] = 668.2143
        ind = ((salt >= 33.915) & (salt < 35))
        mm[ind] = 246.2214
        bb[ind] = -6034.6841
        alkalinity = mm*salt + bb
        return alkalinity
    elif vn == 'phytoplankton':
        phytoplankton = 0.01 + 0*salt
        return phytoplankton
    elif vn == 'chlorophyll':
        chlorophyll = 0.025 + 0*salt
        return chlorophyll
    elif vn == 'zooplankton':
        zooplankton = 0.01 + 0*salt
        return zooplankton
    elif vn == 'SdetritusN':
        SdetritusN = 0 * salt
        return SdetritusN
    elif vn == 'LdetritusN':
        LdetritusN = 0 * salt
        return LdetritusN
    elif vn == 'SdetritusC':
        SdetritusC = 0 * salt
        return SdetritusC
    elif vn == 'LdetritusC':
        LdetritusC = 0 * salt
        return LdetritusC
    
