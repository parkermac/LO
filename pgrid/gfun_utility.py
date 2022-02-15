"""
Utility functions for pgrid.
"""
import numpy as np
import xarray as xr
import seawater as sw

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

def simple_grid(aa, res):
    dlat = aa[3] - aa[2]
    dlon = aa[1] - aa[0]
    mean_lat = np.mean(aa[2:])
    earth_rad = zfun.earth_rad(mean_lat)
    dlatm = earth_rad * np.pi * dlat/180
    dlonm = earth_rad * np.cos(np.pi*mean_lat/180) * np.pi * dlon/180
    nx = int(np.ceil(dlonm/res))
    ny = int(np.ceil(dlatm/res))
    Lon_vec = np.linspace(aa[0], aa[1], nx)
    Lat_vec = np.linspace(aa[2], aa[3], ny)
    return Lon_vec, Lat_vec

def stretched_grid(lon_list, x_res_list, lat_list, y_res_list):
    """
    Input:
    The four lists (units = degrees and meters) are points that define
    segmented vectors along which the resolution changes linearly.
    
    Output:
    vectors of lon and lat suitable for using with meshgrid to make
    lon_rho and lat_rho.
    """
    Lon_list = []
    Lat_list = []
    if (len(lon_list) != len(x_res_list)) or (len(lat_list) != len(y_res_list)):
        print('Lists must be the same length')
        return np.array(Lon_list), np.array(Lat_list)

    lon_vec = np.array(lon_list)
    x_res_vec = np.array(x_res_list)
    lat_vec = np.array(lat_list)
    y_res_vec = np.array(y_res_list)

    R = zfun.earth_rad(np.mean(lat_vec))
    clat = np.cos(np.pi*np.mean(lat_vec)/180)

    lon = lon_list[0]
    Lon_list.append(lon)
    while lon <= lon_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lon]), lon_vec)
        xres = (1-fr)*x_res_vec[i0] + fr*x_res_vec[i1]
        dlon = 180 * xres / (clat * R * np.pi)
        lon = lon + dlon
        Lon_list.append(lon[0])
    lat = lat_list[0]
    Lat_list.append(lat)
    while lat <= lat_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lat]), lat_vec)
        yres = (1-fr)*y_res_vec[i0] + fr*y_res_vec[i1]
        dlat = 180 * yres / (R * np.pi)
        lat = lat + dlat
        Lat_list.append(lat[0])
    return np.array(Lon_list), np.array(Lat_list)

def load_bathy_nc(t_fn):
    ds = xr.open_dataset(t_fn)
    tlon_vec = ds['lon'].values
    tlat_vec = ds['lat'].values
    tz = ds['z'].values
    # There is a bug in xarray with these files: it does
    # not set masked regions to nan.  So we do it by hand.
    tz[tz>1e6] = np.nan
    ds.close()
    return tlon_vec, tlat_vec, tz
    
def load_bathy2(t_fn, lon_vec, lat_vec):
    # load a section of the new NetCDF Smith-Sandwell bathy
    ds = xr.open_dataset(t_fn)
    Lon = ds['lon'].values
    Lat = ds['lat'].values
    i0 = zfun.find_nearest_ind(Lon, lon_vec[0])
    i1 = zfun.find_nearest_ind(Lon, lon_vec[-1])
    j0 = zfun.find_nearest_ind(Lat, lat_vec[0])
    j1 = zfun.find_nearest_ind(Lat, lat_vec[-1])
    tlon_vec = Lon[i0-2:i1+3]
    tlat_vec = Lat[j0-2:j1+3]
    tz = ds['z'][j0-2:j1+3, i0-2:i1+3]
    ds.close()
    return tlon_vec, tlat_vec, tz

def make_nc(out_fn, lon, lat, z, dch):
    """
    Initial creation of the NetCDF grid file.
    """
    
    # create dx, dy (used to make pm and pn)
    plon, plat = pfun.get_plon_plat(lon, lat)
    ulon = plon[:-1, :]
    vlat = plat[:, :-1]
    R = zfun.earth_rad(np.mean(plat[:,1]))
    dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(ulon, axis=1)/180)
    dy = R * (np.pi*np.diff(vlat, axis=0)/180)
    
    # populate dicts of data fields, organized by their dimensions
    #
    # these are on the rho grid
    rho_dict = {
        'h': -z,
        'pm': 1/dx,
        'pn': 1/dy,
        'f': sw.f(lat),
        }
    # these are scalars or characters
    misc_dict = {
        'xl': dx[0,:].sum(),
        'el': dy[:,0].sum(),
        'spherical': 'T',
        }
    # these are on all the grids in tag_list below
    lon_lat_dict = {
        'lon_rho':lon,
        'lat_rho': lat,
        'lon_u': lon[:, :-1] + np.diff(lon, axis=1)/2,
        'lat_u': lat[:, :-1],
        'lon_v': lon[:-1, :],
        'lat_v': lat[:-1, :] + np.diff(lat, axis=0)/2,
        'lon_psi': plon[1:-1, 1:-1],
        'lat_psi': plat[1:-1, 1:-1],
        }
    
    # now populate a data dict with all of these, including dimensions,
    # to feed to an xarray Dataset
    data_dict = dict()
    for vn in rho_dict.keys():
        data_dict[vn] = (('eta_rho', 'xi_rho'), rho_dict[vn])
    for vn in misc_dict.keys():
        data_dict[vn] = (misc_dict[vn])
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        data_dict['lon_'+tag] = (('eta_'+tag, 'xi_'+tag), lon_lat_dict['lon_'+tag])
        data_dict['lat_'+tag] = (('eta_'+tag, 'xi_'+tag), lon_lat_dict['lat_'+tag])
        # masks follow ROMS convention: 1 = water, 0 = land
        data_dict['mask_'+tag] = (('eta_'+tag, 'xi_'+tag), np.ones(lon_lat_dict['lon_'+tag].shape))
    
    # create the Dataset and write to NetCDF (overwrites existing file)
    ds = xr.Dataset(data_dict)
    ds.to_netcdf(out_fn)
    
def make_nudgcoef(dch, out_dir, N, NR, NC):
    # Using info from From https://www.myroms.org/projects/src/ticket/627
    
    out_fn = out_dir / 'nudgcoef.nc'

    fld2 = np.zeros((NR, NC))
    nn = 6 # size of nudging band in gridpoints
    days_short = dch['nudging_days'][0]
    days_long = dch['nudging_days'][1]
    # make inverse time scales (d-1)
    t0 = 1/days_short
    t1 = 1/days_long
    
    if 'north' in dch['nudging_edges']:            
        for i in range(NC):
            for j in range(NR - nn, NR):
                jj = j - NR + nn + 1
                tnud = t1 + jj*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'south' in dch['nudging_edges']:            
        for i in range(NC):
            for j in range(nn):
                tnud = t0 - j*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'west' in dch['nudging_edges']:
        for i in range(nn):
            for j in range(NR):
                tnud = t0 - i*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'east' in dch['nudging_edges']:
        for i in range(NC - nn, NC):
            for j in range(NR):
                ii = i - NC + nn + 1
                tnud = t1 + ii*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])

    fld3 = np.ones((N, 1, 1)) * fld2.reshape((1, NR, NC))
    
    ds = xr.Dataset()
    ds['M2_NudgeCoef'] = (('eta_rho', 'xi_rho'), fld2)
    ds['M2_NudgeCoef'].attrs['long_name'] = '2D momentum inverse nudging coefficients'
    ds['M2_NudgeCoef'].attrs['units'] = 'day-1'
    
    vn_dict = {'M3_NudgeCoef': '3D momentum inverse nudging coefficients',
               'tracer_NudgeCoef': 'generic tracer inverse nudging coefficients',
               'temp_NudgeCoef': 'temp inverse nudging coefficients',
               'salt_NudgeCoef': 'salt inverse nudging coefficients'}
    for vn in vn_dict.keys():
        ds[vn] = (('s_rho', 'eta_rho', 'xi_rho'), fld3)
        ds[vn].attrs['long_name'] = vn_dict[vn]
        ds[vn].attrs['units'] = 'day-1'
    
    for vn in ds.data_vars:
        print('%s %s %s' % (vn, str(ds[vn].shape), ds[vn].attrs['long_name']))
    
    encoding_dict = {vn:{"zlib": True, "complevel": 9} for vn in vn_dict.keys()}
    ds.to_netcdf(out_fn, encoding=encoding_dict)
    
    ds.close()

def GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix,
    fjord_cliff_edges = True, shift=0):
    """
    This is a faster version of GRID_PlusMinusScheme_rx0_ORIG, with about 15x
    speedup in the 100x200 test grid.  It is comparable to the Matlab version.

    ** The depth matrix Hobs MUST BE POSITIVE in non-masked cells **

    The results were nearly identical to those from GRID_PlusMinusScheme_rx0,
    but with some variation up to +/- 45 m in some grid cells.  I suspect that
    this is due to the fact that the order in which I flip the grid around is
    different than in the original.  Since I see no reason for this order
    to be important I will assume the difference is not important.
    
    With fjord_cliff_edges=True is deviates from its usual volume-conserving
    nature when it is next to a masked region, and instead adjusts the slope
    by preferentially deepening at the coast.  This does a much better job of
    preserving thalweg depth in channels like Hood Canal.
    
    """
       
    HH=Hobs.copy()
    HH = HH - shift
    AA = AreaMatrix.copy()
    MM = MSK.copy()
    R=(1-rx0max)/(1+rx0max)
    tol=0.000001
    IsFinished = 1
    count = 0
    maxcount = 1000
    while True and count < maxcount:
        IsFinished=1
        for ff in range(5):
            # the "5" appears to refer to various left-right and up-down
            # flipping permutations 
            if ff == 0:
                do_smooth = True
            elif ff == 1:
                do_smooth = True
                HH = np.fliplr(HH)
                AA = np.fliplr(AA)
                MM = np.fliplr(MM)
            elif ff == 2:
                do_smooth = True
                HH = HH.T
                AA = AA.T
                MM = MM.T
            elif ff == 3:
                do_smooth = True
                HH = np.fliplr(HH)
                AA = np.fliplr(AA)
                MM = np.fliplr(MM)
            elif ff == 4:
                do_smooth = False
                HH = HH.T
                HH = np.fliplr(HH)
                HH = np.flipud(HH)
                AA = AA.T
                AA = np.fliplr(AA)
                AA = np.flipud(AA)
                MM = MM.T
                MM = np.fliplr(MM)
                MM = np.flipud(MM)
            if do_smooth:
                NR, NC = HH.shape
                for ii in range(NC-1):
                    H = HH[:, ii]
                    Hn = HH[:, ii+1]
                    A = AA[:, ii]
                    An = AA[:, ii+1]
                    M = MM[:, ii]
                    Mn = MM[:, ii+1]
                    LowerBound = Hn*R
                    # mask is true when Hn is significantly deeper than H
                    # and when both are water points
                    # and when these are the case it makes H a little deeper
                    # and Hn a litte shallower
                    mask = (H - LowerBound < -tol) & (M == 1) & (Mn == 1)
                    if np.any(mask):
                        IsFinished=0
                        h = (R*Hn - H)/(An + R*A)
                        if ii > 0 and fjord_cliff_edges:
                            Mm = MM[:, ii-1]
                            xm = Mm == 0 # true when there is land to the left
                            xH = H.copy()
                            xH[xm] = xH[xm] + 2*An[xm]*h[xm]
                            xH[~xm] = xH[~xm] + An[~xm]*h[~xm]
                            H = xH.copy()
                            xHn = Hn.copy()
                            xHn[xm] = xHn[xm] - 0*A[xm]*h[xm]
                            xHn[~xm] = xHn[~xm] - A[~xm]*h[~xm]
                            Hn = xHn.copy()
                        else:
                            H = H + An*h
                            Hn = Hn - A*h
                        HH[mask, ii] = H[mask]
                        HH[mask, ii + 1] = Hn[mask]
                        
        if IsFinished == 1:
            break
        #print(count)
        count += 1
    print('Number of iterations = ' + str(count))
    if count == maxcount:
        print('\n** WARNING: more iterations needed! **\n')
    HH = HH + shift
    return HH
    
