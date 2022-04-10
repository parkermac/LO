"""
Functions to use with atmospheric forcing.  Translated from matlab to python.
"""

import numpy as np
import netCDF4 as nc
import seawater as sw
import matplotlib.path as mpath

invar_list = ['Q2', 'T2', 'PSFC', 'U10', 'V10','RAINCV', 'RAINNCV', 'SWDOWN', 'GLW']

outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']

# for plotting, paired with outvar_list
lim_list = [(800,1050), (0,1e-6), (0,900), (0,400), (0,20), (50,100), (-10, 10), (-10,10)]

longname_dict = dict()
units_dict = dict()
timename_dict = dict()
for vn in outvar_list:
    if vn == 'Pair':
        nclongname = 'surface air pressure' 
        ncunits = 'millibar' 
        nctimename = 'pair_time' 
    elif vn == 'rain':
        nclongname = 'rain fall rate' 
        ncunits = 'kilograms meter-2 second-1' 
        nctimename = 'rain_time' 
    elif vn == 'swrad':
        nclongname = 'solar shortwave radiation flux' 
        ncunits = 'watts meter-2' 
        nctimename = 'srf_time' 
    elif vn == 'lwrad_down':
        nclongname = 'downwelling longwave radiation flux' 
        ncunits = 'watts meter-2' 
        nctimename = 'lrf_time' 
    elif vn == 'Tair':
        nclongname = 'surface air temperature' 
        ncunits = 'Celsius' 
        nctimename = 'tair_time' 
    elif vn == 'Qair':
        nclongname = 'surface air relative humidity' 
        ncunits = 'percentage' 
        nctimename = 'qair_time' 
    elif vn == 'Uwind':
        nclongname = 'surface u-wind component' 
        ncunits = 'meter second-1' 
        nctimename = 'wind_time' 
    elif vn == 'Vwind':
        nclongname = 'surface v-wind component' 
        ncunits = 'meter second-1' 
        nctimename = 'wind_time'
    longname_dict[vn] = nclongname
    units_dict[vn] = ncunits
    timename_dict[vn] = nctimename
    
def get_wrf_grid(fn):
    wds = nc.Dataset(fn)
    lon = wds['XLONG'][:].squeeze()
    lat = wds['XLAT'][:].squeeze()
    if False: # Make True to see wrf variable names
        print('\n' + fn.split('/')[-1].center(60,'-'))
        vn_list = []
        for vn in wds.variables:
            vn_list.append(vn)
        print(vn_list)
    wds.close()
    # Get grid size info at the middle of each wrf domain
    NR, NC = lon.shape
    jj = int(NR/2); ii = int(NC/2)
    dx_km, dd_deg = sw.dist(lat[jj,ii], [lon[jj,ii], lon[jj+1,ii+1]])
    return lon, lat, dx_km
    
def get_angle(lon, lat):
    # used to find the angle "theta" of the wrf grids and then return
    # cosine and sine arrays that we use later to rotate the wrf
    # velocities to ROMS E-N orientation.
    NR, NC = lon.shape
    theta = np.nan * np.ones_like(lon)
    for jj in range(NR):
        junk, theta[jj,:-1] = sw.dist(lat[jj,:], lon[jj,:])
    # we repeat the last column because sw.dist returns NC-1
    theta[:,-1] = theta[:,-2]
    if False:
        print(' THETA '.center(60,'-'))
        print(theta[10,-10:])
    # We use negative of the angle because re are rotating back to zero.
    ca = np.cos(-np.pi*theta/180)
    sa = np.sin(-np.pi*theta/180)
    return ca, sa
    
def get_indices_in_polygon(plon_poly, plat_poly, lon, lat):
    # get Boolean mask array "M" that is true for points
    # in lon, lat that are in the polygon plon_poly, plat_poly
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    Rlon = lon.flatten()
    Rlat = lat.flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    M = P.contains_points(R) # boolean
    M = M.reshape(lon.shape)
    return M
    
def gather_and_process_fields(fn, imax, ca, sa, outvar_list):
    # This is where we define any transformations to get from WRF to ROMS variables.
    # we pass outvar_list only because it may have been shortened by the calling program
    # while testing.
    ds = nc.Dataset(fn)
    iv_dict = dict()
    for ivn in invar_list:
        # we trim fields to match the trimmed coordinate arrays
        iv_dict[ivn] = ds[ivn][0,:,:imax].squeeze()
    ds.close()
    # then convert to ROMS units/properties, still on the WRF grid
    # invar_list = ['Q2', 'T2', 'PSFC', 'U10', 'V10','RAINCV', 'RAINNCV', 'SWDOWN', 'GLW']
    # outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']
    ov_dict = dict()
    for ovn in outvar_list:
        if ovn == 'Pair':
            # convert Pa to mbar
            ov_dict[ovn] = iv_dict['PSFC']/100 
        elif ovn == 'rain':
            # set this to zero because (a) we don't really understand the units
            # and (b) is it not used in the simulations at this point 2019.05.22
            ov_dict[ovn] = 0 * (iv_dict['RAINCV']+iv_dict['RAINNCV'])
        elif ovn == 'Tair':
            # convert K to C
            ov_dict[ovn] = iv_dict['T2'] - 273.15
        elif ovn == 'swrad':
            # account for reflection
            ov_dict[ovn] = iv_dict['SWDOWN'] * (1 - 0.1446)
        elif ovn == 'lwrad_down':
            # account for reflection
            ov_dict[ovn] = iv_dict['GLW']
        elif ovn == 'Qair':
            # calculate relative humidity [%]
             ov_dict[ovn] = Z_wmo_RH(ov_dict['Pair'], ov_dict['Tair'], iv_dict['Q2'])
        elif ovn == 'Uwind':
            # % rotate velocity to E-W and N-S
            ov_dict[ovn] = ca*iv_dict['U10'] + sa*iv_dict['V10']
        elif ovn == 'Vwind':
            # % rotate velocity to E-W and N-S
            ov_dict[ovn] = ca*iv_dict['V10'] - sa*iv_dict['U10']
    return ov_dict
    
def interp_to_roms(ov_dict, outvar_list, IMn, NR, NC):
    ovi_dict = dict()
    for ovn in outvar_list:
        v = ov_dict[ovn].flatten()
        ovi_dict[ovn] = v[IMn].reshape((NR,NC))
    return ovi_dict

def Z_wmo_RH(P,T,Q):
    # 5/21/2011 Nick Lederer, modified by Parker MacCready, and recoded
    # from matlab to python by PM 2019.05.16.  Tested against the matlab version
    # using Z_wmo_RH(1000, 10, .001) and both give 13.031628710406915.
    # 
    #  this converts mixing ratio (kg kg-1) which is the usual WRF output [CHECK!], into
    #  relative humidity (%) which is what ROMS expects
    # 
    #  INPUT:
    #  P in hectaPascal or millibar
    #  T in Celcius
    #  Q in kg kg-1
    # 
    #  OUTPUT:
    #  RH in percent
    #
    #  all equations come from Chapter 4 of
    #  http://www.wmo.int/pages/prog/www/IMOP/publications/CIMO-Guide/
    #  CIMO_Guide-7th_Edition-2008.html
    # 
    #  WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION
    #           WMO-No. 8 (Seventh edition) (6 August 2008)
    #
    # Note 2019.05.15: this document no longer exists on the web.

    e_prime = Q*P/(0.62198+Q) # WHO equation 4.A.6
    fp = 1.0016 + 3.15e-6*P - 0.074/P # from Annex 4.B
    ew = 6.112*np.exp(17.62*T/(243.12 + T)) # from Annex 4.B
    ew_prime = fp*ew # from Annex 4.B
    RH = 100 * e_prime/ew_prime # from Annex 4.B

    return RH
    
    
if __name__ == '__main__':
    # examples of uses of the functions, executed in ipython as:
    # run atm_fun.py

    RH = Z_wmo_RH(1000, 10, .001)
    RH_expect = 13.031628710406915
    print('\nTest of Z_wmo_RH:')
    print(' Value - Expected Value = %0.1g\n' % (RH_expect - RH))
    
    print('\nTest of variable attributes')
    for vn in outvar_list:
        print('\n' + vn + ':')
        print(' longname = %s' % (longname_dict[vn]))
        print(' units = %s' % (units_dict[vn]))
        print(' timename = %s' % (timename_dict[vn]))
        