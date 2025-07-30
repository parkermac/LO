"""
Functions to use with atmospheric forcing.  Translated from matlab to python.
"""

import numpy as np
import xarray as xr
import seawater as sw
import matplotlib.path as mpath

invar_list = ['Q2', 'T2', 'PSFC', 'U10', 'V10','RAINC', 'RAINNC', 'SWDOWN', 'GLW']

outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']

    
def get_wrf_grid(fn):
    wds = xr.open_dataset(fn)
    lon = wds['XLONG'].values.squeeze()
    lat = wds['XLAT'].values.squeeze()
    if False: # Make True to see wrf variable names
        print('\n' + fn.name.split('/')[-1].center(60,'-'))
        vn_list = []
        for vn in wds.data_vars:
            print(vn)
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
    ds = xr.open_dataset(fn)
    iv_dict = dict()
    for ivn in invar_list:
        # we trim fields to match the trimmed coordinate arrays
        iv_dict[ivn] = ds[ivn][0,:,:imax].values.squeeze()
    ds.close()
    # then convert to ROMS units/properties, still on the WRF grid
    ov_dict = dict()
    for ovn in outvar_list:
        if ovn == 'Pair':
            # convert Pa to mbar
            ov_dict[ovn] = iv_dict['PSFC']/100 
        elif ovn == 'rain':
            # This is accumulated precipitation [mm] and will be converted to
            # the units expected by ROMS [kg m-2 s-1] in the main code after
            # all hours have been gathered into one array (because we have to
            # take a time derivative).
            ov_dict[ovn] = iv_dict['RAINC']+iv_dict['RAINNC']
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
