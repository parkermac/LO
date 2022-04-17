"""
Functions for tpxo extraction and processing.
"""
from lo_tools import pyTMD_functions as tmd_fun
import xarray as xr
import numpy as np
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
from datetime import datetime, timedelta
import sys

def get_tpxo_clip(Ldir, con, time_dt, domain_tup):
    
    """
    This gets the tpxo9 fields for a given region and day,
    including nodal corrections.  These assume that the time reference
    for the tidal phase is 1/1/1992 (TPXO standard) but the nodal corrections
    are for the day in time_dt.
    
    NOTE: currently this code does not handle cases with positive longitude.
    This problem arises because the tpxo9 arrays are packed with lon 0:360
    whereas in the LO grids we use -180:180.  To fix this we would just need the
    initial extraction from tpxo9 to be clever enough to handle crossing 360
    cleanly.
    
    Input:
    Ldir = the usual dict of LO paths
    con = the tidal constituent to work on, e.g. 'm2'
    time_dt = a datetime for the run day to work on, e.g. datetime(2019,7,4)
    domain_tup = (lon0, lon1, lat0, lat1) a tuple of the limits of the extraction area
    
    Output:
    om, lon, lat, plon, plat, h, amp, phase, umajor, uminor, uincl, uphase
    om = frequency of the constituent [rad s-1]
    lon, lat = plaid arrays of extracted fields
    plon, plat = extended versions of lon, lat, used for pcolormesh plotting
    amp = elevation amplitude [m]
    phase = elevation phase [deg, 0-360]
    umajor = current ellipse major axis [m s-1]
    uminor = current ellipse minor axis [m s-1]
    unicl = current ellipse inclination [deg, 0-180]
    uphase = current ellipse phase [deg, 0-360]
    """

    lon0, lon1, lat0, lat1 = domain_tup
    if lon0 > 0 or lon1 > 0:
        print('tpxo_functions.py error: currently positive lon is not supported')
        sys.exit()
    
    # get some grid and masking info
    in_dir = Ldir['data'] / 'tide' / 'tpxo9'
    g_fn = in_dir / 'grid_tpxo9_atlas_30_v5.nc'
    g_ds = xr.open_dataset(g_fn)
    # Get grid depth, lon, and lat for a clipped region,
    # along with clipping indices.
    lon_vec = g_ds.lon_z.values # 0:360
    lat_vec = g_ds.lat_z.values # -90:90
    i0 = zfun.find_nearest_ind(lon_vec, lon0 + 360)
    i1 = zfun.find_nearest_ind(lon_vec, lon1 + 360)
    j0 = zfun.find_nearest_ind(lat_vec, lat0)
    j1 = zfun.find_nearest_ind(lat_vec, lat1)
    lon, lat = np.meshgrid(lon_vec[i0:i1], lat_vec[j0:j1])
    lon = lon - 360 # hack to convert to -180:180 that only works west of 0 degrees
    # extended lon, lat for pcolormesh plotting
    plon, plat = pfun.get_plon_plat(lon, lat)
    # depth field (0 = land)
    h = g_ds.hz[i0:i1, j0:j1].values
    h = h.T

    # load tpxo9 fields
    c_fn = in_dir / ('h_' + con + '_tpxo9_atlas_30_v5.nc')
    c_ds = xr.open_dataset(c_fn)

    u_fn = in_dir / ('u_' + con + '_tpxo9_atlas_30_v5.nc')
    u_ds = xr.open_dataset(u_fn)

    # Load contituent fields
    # Info from c_ds.hRe.field (MATLAB notation):
    #    amp=abs(hRe+i*hIm)
    #    GMT phase=atan2(-hIm,hRe)/pi*180 [converting radians to degrees]
    # Also 0 = land so we use that as a mask.
    # The amplitude units are mm (int), and I believe the phase is relative to 1/1/1992
    hRe = c_ds.hRe[i0:i1, j0:j1].values
    hIm = c_ds.hIm[i0:i1, j0:j1].values
    # Note that these are packed [lon, lat] hence the transpose below to get them
    # to my standard [lat, lon] packing.
    # Real part
    hRe = np.float64(hRe)
    hRe = hRe.T
    # Imaginary part
    hIm = np.float64(hIm)
    hIm = hIm.T
    # Complex
    hcx = hRe + 1j*hIm
    # Then find the amplitude and phase of the complex amplitide
    amp = abs(hcx)
    phase = np.arctan2(-hIm,hRe) # -pi:pi
    phase[h==0] = np.nan
    amp[h==0] = np.nan

    # then get the velocities
    uRe = u_ds.uRe[i0:i1+1, j0:j1].values
    uIm = u_ds.uIm[i0:i1+1, j0:j1].values
    vRe = u_ds.vRe[i0:i1, j0:j1+1].values
    vIm = u_ds.vIm[i0:i1, j0:j1+1].values

    uRe = np.float64(uRe)
    uRe = uRe.T
    uIm = np.float64(uIm)
    uIm = uIm.T
    ucx = uRe + 1j*uIm
    # interpolate to the z grid and mask
    ucx = (ucx[:,:-1] + ucx[:,1:])/2
    # convert from cm2/s to m/s
    ucx[h>0] = (ucx[h>0]/1e4) / h[h>0]
    ucx[h==0] = np.nan

    vRe = np.float64(vRe)
    vRe = vRe.T
    vIm = np.float64(vIm)
    vIm = vIm.T
    vcx = vRe + 1j*vIm
    # interpolate to the z grid and mask
    vcx = (vcx[:-1,:] + vcx[1:,:])/2
    # convert from cm2/s to m/s
    vcx[h>0] = (vcx[h>0]/1e4) / h[h>0]
    vcx[h==0] = np.nan

    # Calculate current ellipse parametters
    (umajor, uminor, uincl, uphase) = tmd_fun.tidal_ellipse(ucx, vcx)
    # note that:
    # uincl is degresse 0:180
    # uphase is segrees 0:360

    # Apply nodal corrections and Greenwich phase using pyTMD functions
    junk, ph, om, junk, junk = tmd_fun.load_constituent(con)
    # ph = Greenwich phase of this constituent [rad] without nodal correction
    # om = frequency of this constituent [rad s-1]
    # 2*np.pi/(om*3600) gives period in hours.
    # The pyTMD nodal correction code expects time in Modified Julian Date (mjd)
    # which is the number of days since midnight on November 17, 1858.  Weird.
    mjd = (time_dt - datetime(1858,11,17)).days
    pu, pf, G = tmd_fun.load_nodal_corrections(mjd, [con])
    # pu = nodal correction of phase [rad]
    # pf = nodal correction to multiply amplitude by [dimensionless]
    # include Greenwich phase and nodal phase corrrection
    phase = phase - ph - pu
    uphase = np.pi*uphase/180 - ph - pu
    # wrap phase to ensure it is in -pi to pi
    phase = np.angle(np.exp(1j*phase))
    uphase = np.angle(np.exp(1j*uphase))
    # convert to degrees [-180:180]
    phase = 180 * phase / np.pi
    uphase = 180 * uphase / np.pi
    # convert to 0-360 format, just so we conform to some convention
    phase[phase<0] += 360 # [0:360]
    uphase[uphase<0] += 360 # [0:360]
    # apply nodal correction to amplitude
    amp = pf * amp / 1000 # also convert mm to m
    umajor = pf * umajor
    uminor = pf * uminor
    
    return om, lon, lat, plon, plat, h, amp, phase, umajor, uminor, uincl, uphase