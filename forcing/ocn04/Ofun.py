"""
Functions for getting and processing the ocean forcing.

"""
import os, sys
from datetime import datetime, timedelta
import numpy as np
import pickle
from scipy.spatial import cKDTree
import gsw
import warnings

import xarray as xr
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import pandas as pd

from lo_tools import Lfun, zfun, zrfun
from lo_tools import hycom_functions as hfun

def messages(mess_str, stdout, stderr):
    # utility function to help with subprocess errors
    try:
        if len(stdout) > 0:
            print(mess_str)
            print(stdout.decode())
    except TypeError:
        pass
    try:
        if len(stderr) > 0:
            print(mess_str)
            print(stderr.decode())
    except TypeError:
        pass

def get_indices(h_out_dir, dt0, dt1, url_dict, verbose=False):
    # find and check the indices into hycom for the extraction

    # specify the sub region of hycom to extract
    aa = hfun.aa
    # convert to hycom format
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360

    ind_dicts = dict()
    for hkey in url_dict.keys():
        got_indices = True
        if verbose:
            print('-- getting indices for: ' + hkey)
            sys.stdout.flush()
        url = url_dict[hkey]
        out_fn = h_out_dir / (hkey + '_tyx.nc')
        # get rid of the old version, if it exists
        out_fn.unlink(missing_ok=True)
        # extract coordinates
        cmd_list = ['ncks','-O','-v','time,lat,lon',url,str(out_fn)]
        #print(cmd_list)
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        messages('get_indices() messages:', stdout, stderr)
        if len(stderr) > 0:
            got_indices = False
            break
        # use the results
        ds = xr.open_dataset(out_fn)
        # find selected indices to use with ncks to extract fields
        t = ds.time.values
        tind = pd.DatetimeIndex(t)
        it0 = np.argwhere(tind==dt0)[0][0]
        it1 = np.argwhere(tind==dt1)[0][0]
        x = ds.lon.values
        y = ds.lat.values
        ix0 = zfun.find_nearest_ind(x,west)
        ix1 = zfun.find_nearest_ind(x,east)
        iy0 = zfun.find_nearest_ind(y,south)
        iy1 = zfun.find_nearest_ind(y,north)
        ds.close()
        ind_dict = {'it0':it0, 'it1':it1, 'ix0':ix0, 'ix1':ix1, 'iy0':iy0, 'iy1':iy1}
        ind_dicts[hkey] = ind_dict

    return ind_dicts, got_indices

def get_hycom_data(out_fn, hkey, ind_dicts, url_dict, hycom_var_dict, verbose=False): # get hourly or 3-hourly data
    """"
    Gets arrays of hycom extractions over some time range.
    """
    # get rid of the old version, if it exists
    out_fn.unlink(missing_ok=True)
    print(' - getting hycom fields for ' + str(out_fn))
    got_hycom_data = True

    tt0 = time()
    ind_dict = ind_dicts[hkey]
    it0 = ind_dict['it0']
    it1 = ind_dict['it1']
    ix0 = ind_dict['ix0']
    ix1 = ind_dict['ix1']
    iy0 = ind_dict['iy0']
    iy1 = ind_dict['iy1']
    url = url_dict[hkey]
    hvar = hycom_var_dict[hkey]
    # extract data from HYCOM file
    cmd_list = ['ncks','-O','-v',hvar,
        '-d','time,'+str(it0)+','+str(it1),
        '-d','lat,'+str(iy0)+','+str(iy1),
        '-d','lon,'+str(ix0)+','+str(ix1),
        url,str(out_fn)]
    if verbose:
        print(cmd_list)
        sys.stdout.flush()
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    messages('ncks extract messages:', stdout, stderr)
            
    if len(stderr) > 0:
        got_hycom_data = False
        return got_hycom_data
    else:
        print('Took %0.1f sec to get %s' % (time()-tt0, hvar))
        return got_hycom_data
            
def get_coords(in_dir):
    """
    get coordinate fields and sizes
    """
    coord_dict = pickle.load(open(in_dir / 'coord_dict.p', 'rb'))
    lon = coord_dict['lon']
    lat = coord_dict['lat']
    z = coord_dict['z']
    L = len(lon)
    M = len(lat)
    N = len(z)
    # Create arrays of distance from the center (m) so that the
    # nearest neighbor extrapolation is based on physical distance
    Lon, Lat = np.meshgrid(lon,lat)
    X, Y = zfun.ll2xy(Lon, Lat, lon.mean(), lat.mean())
    return (lon, lat, z, L, M, N, X, Y)

def checknan(fld):
    """
    A utility function that issues a warning if there are nans in fld.
    """
    if np.isnan(fld).sum() > 0:
        print('WARNING: nans in data field')    

def extrap_nearest_to_masked(X, Y, fld, fld0=0):
    """
    INPUT: fld is a 2D array (np.ndarray or np.ma.MaskedArray) on spatial grid X, Y
    OUTPUT: a numpy array of the same size with no mask
    and no missing values.        
    If input is a masked array:        
        * If it is ALL masked then return an array filled with fld0.         
        * If it is PARTLY masked use nearest neighbor interpolation to
        fill missing values, and then return data.        
        * If it is all unmasked then return the data.    
    If input is not a masked array:        
        * Return the array.    
    """
    # first make sure nans are masked
    if np.ma.is_masked(fld) == False:
        fld = np.ma.masked_where(np.isnan(fld), fld)
        
    if fld.all() is np.ma.masked:
        #print('  filling with ' + str(fld0))
        fldf = fld0 * np.ones(fld.data.shape)
        fldd = fldf.data
        checknan(fldd)
        return fldd
    else:
        # do the extrapolation using nearest neighbor
        fldf = fld.copy() # initialize the "filled" field
        xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
        xynew = np.array((X[fld.mask],Y[fld.mask])).T
        a = cKDTree(xyorig).query(xynew)
        aa = a[1]
        fldf[fld.mask] = fld[~fld.mask][aa]
        fldd = fldf.data
        checknan(fldd)
        return fldd

def get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir):
    """
    Make use of extrap_nearest_to_masked() to fill fields completely
    before interpolating to the ROMS grid.
    
    It also creates ubar and vbar, and converts the temperature to potential temperature.
    """
    b = pickle.load(open(in_fn, 'rb'))
    vn_list = list(b.keys())    
    # check that things are the expected shape
    def check_coords(shape_tuple, arr_shape):
        if arr_shape != shape_tuple:
            print('WARNING: array shape mismatch')
    for vn in vn_list:
        if vn == 'dt':
            pass
        elif vn == 'ssh':
            check_coords((M, L), b[vn].shape)
        else:
            check_coords((N, M, L), b[vn].shape)    
    # creat output array and add dt to it.
    vn_list.remove('dt')
    V = dict()
    for vn in vn_list:
        V[vn] = np.nan + np.ones(b[vn].shape)
    V['dt'] = b['dt']    
    # extrapolate ssh
    vn = 'ssh'
    vssh = b[vn].copy()
    vvssh = extrap_nearest_to_masked(X, Y, vssh)
    V[vn] = vvssh
    vn_list.remove('ssh')    
    # extrapolate 3D fields
    for vn in vn_list:
        v = b[vn].copy()
        if vn == 't3d':
            v0 = np.nanmin(v)
        elif vn == 's3d':
            v0 = np.nanmax(v)
        if vn in ['t3d', 's3d']:
            # print(' -- extrapolating ' + vn)
            for k in range(N):
                fld = v[k, :, :]
                fldf = extrap_nearest_to_masked(X, Y, fld, fld0=v0)
                V[vn][k, :, :] = fldf
        elif vn in ['u3d', 'v3d']:
            # print(' -- extrapolating ' + vn)
            vv = v.copy()
            vv = np.ma.masked_where(np.isnan(vv), vv)
            vv[vv.mask] = 0
            V[vn] = vv.data

    # Create ubar and vbar.
    # Note: this is slightly imperfect because the z levels are at the same
    # position as the velocity levels.
    dz = np.nan * np.ones((N, 1, 1))
    dz[1:, 0, 0]= np.diff(z)
    dz[0, 0, 0] = dz[1, 0, 0]
    
    # account for the fact that the new hycom fields do not show up masked
    u3d = np.ma.masked_where(np.isnan(b['u3d']),b['u3d'])
    v3d = np.ma.masked_where(np.isnan(b['v3d']),b['v3d'])
    dz3 = dz * np.ones_like(u3d) # make dz a masked array
    b['ubar'] = np.sum(u3d*dz3, axis=0) / np.sum(dz3, axis=0)
    b['vbar'] = np.sum(v3d*dz3, axis=0) / np.sum(dz3, axis=0)
    
    for vn in ['ubar', 'vbar']:
        v = b[vn]
        vv = v.copy()
        vv = np.ma.masked_where(np.isnan(vv), vv)
        vv[vv.mask] = 0
        V[vn] = vv.data
        
    # calculate potential temperature
    Lon, Lat = np.meshgrid(lon,lat)
    Z = z.reshape((len(z),1,1))
    p = gsw.p_from_z(Z,Lat)
    SA = gsw.SA_from_SP(V['s3d'],p,Lon,Lat)
    th = gsw.pt0_from_t(SA, V['t3d'], p)
    # debugging
    # for xx in [SA, V['t3d'], p]:
    #     print(np.min(xx))
    #     print(np.max(xx))
    V['theta'] = th

    return V

def get_xyr(G, vn):
    """
    A utility function for getting any of the ROMS grids.
    """
    if vn in ['ssh', 'theta', 's3d']:
        xr = G['lon_rho']
        yr = G['lat_rho']
    elif vn in ['ubar', 'u3d']:
        xr = G['lon_u']
        yr = G['lat_u']
    elif vn in ['vbar', 'v3d']:
        xr = G['lon_v']
        yr = G['lat_v']
    else:
        print('Unknown variable name for get_xyr: ' + vn)
    return xr, yr

def get_zr(G, S, vn):
    """
    A utility function to get the ROMS z coordinate on any of
    the ROMS grids.
    """
    h = G['h']
    if vn in ['theta', 's3d']:
        zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    elif vn in ['u3d']:    
        xru, yru = get_xyr(G, 'ubar')
        hu = (h[:,1:] + h[:,:-1])/2
        zr = zrfun.get_z(hu, 0*hu, S, only_rho=True)    
    elif vn in ['v3d']:    
        hv = (h[1:, :] + h[:-1,:])/2
        zr = zrfun.get_z(hv, 0*hv, S, only_rho=True)
    else:
        print('Unknown variable name for get_zr: ' + vn)
    return zr

def get_zinds(h, S, z):
    """
    Precalculate the array of indices to go from HYCOM z to ROMS z.
    This just finds the vertical index in HYCOM z for each ROMS z_rho
    value in the whole 3D array, with the index being the UPPER one of
    the two HYCOM z indices that any ROMS z falls between.
    """
    tt0 = time()
    zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    zrf = zr.flatten()
    zinds = np.nan * np.ones_like(zrf)
    if isinstance(z, np.ma.MaskedArray):
        z = z.data
    for ii in range(len(z)-1):
        zlo = z[ii]; zhi = z[ii+1]
        mask = (zrf>zlo) & (zrf<=zhi)
        zinds[mask] = ii+1 # this is where the UPPER index is enforced
    zinds = zinds.astype(int)
    if isinstance(zinds, np.ma.MaskedArray):
        zinds = zinds.data
    # if verbose:
    #     print(' --create zinds array took %0.1f seconds' % (time() - tt0))
    return zinds

def get_interpolated(G, S, b, lon, lat, z, N, zinds):
    """
    This does the horizontal and vertical interpolation to get from
    extrapolated, filtered HYCOM fields to ROMS fields.

    We use fast nearest neighbor interpolation as much as possible.
    Also we interpolate everything to the ROMS rho grid, and then crudely
    interpolate to the u and v grids at the last moment.  Much simpler.
    """
    
    # start input dict
    c = {}
    
    # precalculate useful arrays that are used for horizontal interpolation
    if isinstance(lon, np.ma.MaskedArray):
        lon = lon.data
    if isinstance(lat, np.ma.MaskedArray):
        lat = lat.data
    Lon, Lat = np.meshgrid(lon,lat)
    XYin = np.array((Lon.flatten(), Lat.flatten())).T
    XYr = np.array((G['lon_rho'].flatten(), G['lat_rho'].flatten())).T
    h = G['h']
    IMr = cKDTree(XYin).query(XYr)[1]
    
    # 2D fields
    for vn in ['ssh', 'ubar', 'vbar']:
        vv = b[vn].flatten()[IMr].reshape(h.shape)
        if vn == 'ubar':
            vv = (vv[:,:-1] + vv[:,1:])/2
        elif vn == 'vbar':
            vv = (vv[:-1,:] + vv[1:,:])/2
        vvc = vv.copy()
        # always a good idea to make sure dict entries are not just pointers
        # to arrays that might be changed later, hence the .copy()
        c[vn] = vvc
        checknan(vvc)
        
    # 3D fields
    # create intermediate arrays which are on the ROMS lon_rho, lat_rho grid
    # but have the HYCOM vertical grid (N layers)
    F = np.nan * np.ones(((N,) + h.shape))
    vi_dict = {}
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        FF = F.copy()
        for nn in range(N):
            vin = b[vn][nn,:,:].flatten()
            FF[nn,:,:] = vin[IMr].reshape(h.shape)
        checknan(FF)
        vi_dict[vn] = FF
    
    # do the vertical interpolation from HYCOM to ROMS z positions
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        vi = vi_dict[vn]
        hinds = np.indices((S['N'], G['M'], G['L']))
        vvf = vi[zinds, hinds[1].flatten(), hinds[2].flatten()]
        vv = vvf.reshape((S['N'], G['M'], G['L']))
        vvc = vv.copy()
        if vn == 'u3d':
            vvc = (vvc[:,:,:-1] + vvc[:,:,1:])/2
        elif vn == 'v3d':
            vvc = (vvc[:,:-1,:] + vvc[:,1:,:])/2
        checknan(vvc)
        c[vn] = vvc
    return c


