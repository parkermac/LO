"""
Functions for extracting and processing GLORYS ocean fields.
"""

import copernicusmarine
from datetime import datetime, timedelta
from time import time
import xarray as xr
import numpy as np

from lo_tools import Lfun, zfun, zrfun
Ldir = Lfun.Lstart()

# get Copernicus login credentials
cred_fn = Ldir['data'] / 'accounts' / 'glorys_pm_2025.04.02.csv'
cred_dict = Lfun.csv_to_dict(cred_fn)

# set lat-lon limits for the glorys extraction
aa = [-131, -121, 41, 53]

def get_glorys_forecast(dt, out_dir, verbose=False):
    """
    Extract glorys daily average from forecast at one time.
    dt = datetime for the extraction: should be at midnight
    out_dir = Path to output directory
    """
    # NOTE: For the forecast you have to get each variable separately.
    # set variables to get
    vn_list = ['so', 'thetao', 'zos','cur']
    tt00 = time()
    g_dstr = dt.strftime('%Y-%m-%dT%H:%M:%S')
    for vn in vn_list:
        tt0 = time()
        if vn == 'zos':
            dataset_id='cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
        else:
            dataset_id='cmems_mod_glo_phy-'+vn+'_anfc_0.083deg_P1D-m'
        if vn == 'cur':
            variables=['uo','vo']
        else:
            variables=[vn]
        out_name = 'forecast_'+vn+'.nc'
        copernicusmarine.subset(
            dataset_id=dataset_id,
            variables=variables,
            minimum_longitude=aa[0],
            maximum_longitude=aa[1],
            minimum_latitude=aa[2],
            maximum_latitude=aa[3],
            start_datetime=g_dstr,
            end_datetime=g_dstr,
            output_filename=out_name,
            output_directory=str(out_dir),
            username=cred_dict['username'],
            password=cred_dict['password'],
            overwrite=True,
            disable_progress_bar=True
        )
        if verbose:
            print('\n - Time to get %s at one time = %0.1f' % (vn, time()-tt0))
        if verbose:
            out_fn = out_dir / out_name
            ds = xr.open_dataset(out_fn)
            print(ds)
            ds.close()
    if verbose:
        print('\nTime to get all variables from forecast = %0.1f' % (time()-tt00))

def get_glorys_hindcast(dt, out_dir, verbose=False):
    """
    Extract glorys daily average from hindcast at one time.
    dt = datetime for the extraction: should be at midnight
    out_dir = Path to output directory
    """
    # NOTE: For the hindcast you can get all variables at once.
    # set variables to get
    vn_list = ['uo', 'vo', 'so', 'thetao', 'zos']
    dstr = dt.strftime(Lfun.ds_fmt)
    out_name = 'hindcast_' + dstr + '.nc'
    g_dstr = dt.strftime('%Y-%m-%dT%H:%M:%S')
    tt0 = time()
    copernicusmarine.subset(
        dataset_id='cmems_mod_glo_phy_my_0.083deg_P1D-m',
        variables=vn_list,
        minimum_longitude=aa[0],
        maximum_longitude=aa[1],
        minimum_latitude=aa[2],
        maximum_latitude=aa[3],
        start_datetime=g_dstr,
        end_datetime=g_dstr,
        output_filename=out_name,
        output_directory=str(out_dir),
        username=cred_dict['username'],
        password=cred_dict['password'],
        overwrite=True,
        disable_progress_bar=True
    )
    if verbose:
        print('\nTime to get all variables from hindcast = %0.1f' % (time()-tt0))

def create_roms_grid_info(gridname):
    """
    Create grid info for a specific grid, e.g. gridname = 'cas7'
    """
    grid_dir = Ldir['data'] / 'grids' / gridname
    G = zrfun.get_basic_info(grid_dir / 'grid.nc', only_G=True)
    S_info_dict = Lfun.csv_to_dict(grid_dir / 'S_COORDINATE_INFO.csv')
    S = zrfun.get_S(S_info_dict)
    return G, S

def get_zr(G, S, vn):
    """
    A utility function to get the ROMS z coordinate on any of
    the ROMS grids.
    """
    h = G['h']
    if vn in ['salt', 'temp', 'zeta']:
        zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    elif vn == 'u':    
        hu = (h[:,1:] + h[:,:-1])/2
        zr = zrfun.get_z(hu, 0*hu, S, only_rho=True)    
    elif vn == 'v':    
        hv = (h[1:, :] + h[:-1,:])/2
        zr = zrfun.get_z(hv, 0*hv, S, only_rho=True)
    else:
        print('Unknown variable name for get_zr: ' + vn)
    return zr

def interpolate_glorys_to_roms(fng, vn, vng, gtag, zr, G, verbose=False):
    """
    This function interpolates a glorys field to a roms grid. First it
    fills what it can using linear interpolation from the regular glorys
    grid (fast), and then using cKDTree nearest-neighbor to fill in any
    remaining gaps (slower, about 8 seconds per 3-D field).

    This function is specific to 3-D fields.
    """
    # open the raw glorys field
    dsg = xr.open_dataset(fng)
    tt0 = time()
    # Interpolate glorys field to our ROMS grid.
    # start by using linear interpolation...
    from scipy.interpolate import RegularGridInterpolator
    xx = dsg.longitude.to_numpy()
    yy = dsg.latitude.to_numpy()
    zz = - dsg.depth.to_numpy()
    data = dsg[vng][0,:,:,:].to_numpy()
    interp = RegularGridInterpolator((zz,yy,xx), data,
        method='linear', bounds_error=False)
    # points in the ROMS grid
    N,M,L = zr.shape
    X = np.tile(G['lon_'+gtag].reshape(1,M,L),[N,1,1])
    Y = np.tile(G['lat_'+gtag].reshape(1,M,L),[N,1,1])
    Z = zr
    mask = G['mask_'+gtag]==1
    Mask = np.tile(mask.reshape(1,M,L),[N,1,1])
    zyx = np.array((Z[Mask].flatten(),Y[Mask].flatten(),X[Mask].flatten())).T
    interpolated_values = interp(zyx)
    FLD = np.nan * zr
    FLD[Mask] = interpolated_values
    if verbose:
        print('- time to interpolate = %0.1f sec' % (time()-tt0))
    if True:
        # since this is the slowest step you can turn it off for testing
        tt0 = time()
        # Next fill in remaining missing values using nearest neighbor.
        from scipy.spatial import cKDTree
        z2,y2,x2 = np.meshgrid(zz,yy,xx, indexing='ij')
        mask2 = ~ np.isnan(data)
        Data = data[mask2].flatten()
        zyx2 = np.array((z2[mask2].flatten(),y2[mask2].flatten(),x2[mask2].flatten())).T
        zyxT = cKDTree(zyx2)
        print('- time to make tree = %0.1f sec' % (time()-tt0))
        tt0 = time()
        mask3 = np.isnan(FLD) & Mask
        zyx3 = np.array((Z[mask3].flatten(),Y[mask3].flatten(),X[mask3].flatten())).T
        fill_data = Data[zyxT.query(zyx3, workers=-1)[1]]
        FLD[mask3] = fill_data
        if verbose:
            print('- time to use tree = %0.1f sec' % (time()-tt0))
        # RESULT: Using the tree is the slowest part, 8 sec per 3-D field on my mac.
        # It can slow down in a long interactive ipython session, maybe a garbage
        # collection issue?
    dsg.close()
    return FLD

def interpolate_glorys_to_roms_2d(fng, vn, vng, gtag, G, verbose=False):
    """
    This function interpolates a glorys field to a roms grid. First it
    fills what it can using linear interpolation from the regular glorys
    grid (fast), and then using cKDTree nearest-neighbor to fill in any
    remaining gaps (slower, about 8 seconds per 3-D field).

    This function is specific to 2-D fields, presumably just z0s => zeta.
    """
    # open the raw glorys field
    dsg = xr.open_dataset(fng)
    tt0 = time()
    # Interpolate glorys field to our ROMS grid.
    # start by using linear interpolation...
    from scipy.interpolate import RegularGridInterpolator
    xx = dsg.longitude.to_numpy()
    yy = dsg.latitude.to_numpy()
    data = dsg[vng][0,:,:].to_numpy()
    interp = RegularGridInterpolator((yy,xx), data,
        method='linear', bounds_error=False)
    # points in the ROMS grid
    X = G['lon_'+gtag]
    Y = G['lat_'+gtag]
    Mask = G['mask_'+gtag]==1
    yx = np.array((Y[Mask].flatten(),X[Mask].flatten())).T
    interpolated_values = interp(yx)
    FLD = np.nan * G['h']
    FLD[Mask] = interpolated_values
    if verbose:
        print('- time to interpolate = %0.1f sec' % (time()-tt0))
    if True:
        # since this is the slowest step you can turn it off for testing
        tt0 = time()
        # Next fill in remaining missing values using nearest neighbor.
        from scipy.spatial import cKDTree
        y2,x2 = np.meshgrid(yy,xx, indexing='ij')
        mask2 = ~ np.isnan(data)
        Data = data[mask2].flatten()
        yx2 = np.array((y2[mask2].flatten(),x2[mask2].flatten())).T
        yxT = cKDTree(yx2)
        print('- time to make tree = %0.1f sec' % (time()-tt0))
        tt0 = time()
        mask3 = np.isnan(FLD) & Mask
        yx3 = np.array((Y[mask3].flatten(),X[mask3].flatten())).T
        fill_data = Data[yxT.query(yx3, workers=-1)[1]]
        FLD[mask3] = fill_data
        if verbose:
            print('- time to use tree = %0.1f sec' % (time()-tt0))
        # RESULT: Very fast.
    dsg.close()
    return FLD
