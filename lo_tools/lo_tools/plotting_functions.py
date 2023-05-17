"""
Module of basic utilities for plotting.  Used heavily in plotting/roms_plots.py
"""

import sys
from pathlib import Path
import pandas as pd

# get initial version of Ldir when this module is loaded
from lo_tools import Lfun, zrfun, zfun
Ldir = Lfun.Lstart()

import numpy as np
from datetime import datetime, timedelta
import pytz

if ('_mac' in Ldir['lo_env']) or ('_pc' in Ldir['lo_env']): # mac or pc version
    #print('mac version') # debugging
    pass
else: # remote/linux version
    #print('linux version')
    import matplotlib as mpl
    mpl.use('Agg')
# do these after the mpl.use() call to get remote plotting to work
import matplotlib.pyplot as plt
import matplotlib.path as mpath

# some handy dictionaries
month_name_dict = dict(zip(range(1,13), ['January', 'February', 'March',
    'April', 'May', 'June', 'July', 'August', 'September', 'October',
    'November', 'December']))
month_color_dict = dict(zip(range(1,13),
    ['mediumblue', 'royalblue', 'cadetblue', 'aquamarine',
    'lightgreen', 'greenyellow', 'gold', 'orange',
    'lightsalmon', 'mediumorchid', 'slateblue', 'purple']))
bbox = dict(facecolor='w', edgecolor='None',alpha=.5)

def get_bot_top_arag(fn, aa=[]):
    """
    Returns bottom and top Aragonite Saturation State fields for a given history file.
    
    You can also provide an optional list of axis limits, "aa".
    
    It also returns plon and plat for pcolormesh plotting.
    """
    import xarray as xr
    from PyCO2SYS import CO2SYS
    import gsw
    G = zrfun.get_basic_info(fn, only_G=True)
    ds = xr.open_dataset(fn)
    if len(aa) == 4:
        # find indices that encompass region aa
        i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0])
        i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1])
        j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2])
        j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3])
        # print('%d %d %d %d' % (i0, i1, j0, j1))
        # print(str(G['lon_rho'].shape))
        # print(str(G['lon_rho'][j0:j1,i0:i1].shape))
    else:
        nrows, ncols = G['lon_rho'].shape
        i0 = 0; i1=ncols
        j0 = 0; j1 = nrows
    lon = G['lon_rho'][j0:j1,i0:i1]
    lat = G['lat_rho'][j0:j1,i0:i1]
    hh = G['h'][j0:j1,i0:i1]
    px, py = get_plon_plat(lon, lat)
    v_dict = dict()
    vn_in_list = ['temp', 'salt', 'alkalinity', 'TIC']
    for nlev in [0,-1]:
        if nlev == 0:
            zz = -hh
        elif nlev == -1:
            zz = 0 * hh
        for cvn in vn_in_list:
            v_dict[cvn] = ds[cvn][0,nlev,j0:j1,i0:i1].values.squeeze()
            # # debugging
            # print('%s %s' % (cvn, str(v_dict[cvn].shape)))
            # sys.stdout.flush()
        pres = gsw.p_from_z(zz, lat) # pressure [dbar]
        SA = gsw.SA_from_SP(v_dict['salt'], pres, lon, lat)
        CT = gsw.CT_from_pt(SA, v_dict['temp'])
        rho = gsw.rho(SA, CT, pres) # in situ density
        temp = gsw.t_from_CT(SA, CT, pres) # in situ temperature
        # convert from umol/L to umol/kg using in situ dentity
        alkalinity = 1000 * v_dict['alkalinity'] / rho
        alkalinity[alkalinity < 100] = np.nan
        TIC = 1000 * v_dict['TIC'] / rho
        TIC[TIC < 100] = np.nan
        # See LPM/co2sys_test/test0.py for info.
        import PyCO2SYS as pyco2
        CO2dict = pyco2.sys(par1=alkalinity, par2=TIC, par1_type=1, par2_type=2,
            salinity=v_dict['salt'], temperature=temp, pressure=pres,
            total_silicate=50, total_phosphate=2,
            opt_pH_scale=1, opt_k_carbonic=10, opt_k_bisulfate=1)
        ph = CO2dict['pH']
        arag = CO2dict['saturation_aragonite']
        ph = ph.reshape(pres.shape)
        arag = arag.reshape(pres.shape)
        # account for WET_DRY
        if 'wetdry_mask_rho' in ds.data_vars:
            mwd = ds.wetdry_mask_rho[0,j0:j1,i0:i1].values.squeeze()
            arag[mwd==0] = np.nan
        if nlev == 0:
            arag_bot = arag
        elif nlev == -1:
            arag_top = arag
    return arag_bot, arag_top, px, py

def start_plot(fs=14, figsize=(14,10)):
    plt.rc('font', size=fs)
    plt.rc('figure', figsize=figsize)
    
def end_plot():
    plt.rcdefaults()

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.cos(np.pi*yav/180))

def add_coast(ax, dir0=Ldir['data'], color='k', linewidth=0.5):
    fn = str(dir0) + '/coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    ax.plot(C['lon'].values, C['lat'].values, '-', color=color, linewidth=linewidth)

def get_coast(dir0=Ldir['data']):
    fn = str(dir0) + '/coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    lon = C['lon'].values
    lat = C['lat'].values
    return lon, lat
    
def get_plon_plat(lon, lat):
    """
    This takes the 2-D lon and lat grids (ndarrays) and returns extended
    "psi" grids that are suitable for plotting using pcolormesh
    for any field on the original grid.
    NOTE: It checks to make sure the original grid is plaid.
    
    You would pass it G['lon_rho'] and G['lat_rho'], for
    a field on the rho grid.
    """
    # input checking
    Lon = lon[0,:]
    Lat = lat[:,0]
    if (Lon - lon[-1,:]).sum() != 0:
        print('Error from get_plon_plat: lon grid not plaid')
        sys.exit()
    if (Lat - lat[:,-1]).sum() != 0:
        print('Error from get_plon_plat: lat grid not plaid')
        sys.exit()
    plon = np.ones(len(Lon) + 1)
    plat = np.ones(len(Lat) + 1)
    dx2 = np.diff(Lon)/2
    dy2 = np.diff(Lat)/2
    Plon = np.concatenate(((Lon[0]-dx2[0]).reshape((1,)), Lon[:-1]+dx2, (Lon[-1]+dx2[-1]).reshape((1,))))
    Plat = np.concatenate(((Lat[0]-dy2[0]).reshape((1,)), Lat[:-1]+dy2, (Lat[-1]+dy2[-1]).reshape((1,))))
    plon, plat = np.meshgrid(Plon, Plat)
    return plon, plat

def auto_lims(fld, vlims_fac=3):
    """
    A convenience function for automatically setting color limits.
    Input: a numpy array (masked OK)
    Output: tuple of good-guess colorscale limits for a pcolormesh plot.    
    """
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    flo = np.nanmax([np.nanmean(fld) - vlims_fac*np.nanstd(fld), np.nanmin(fld)])
    fhi = np.nanmin([np.nanmean(fld) + vlims_fac*np.nanstd(fld), np.nanmax(fld)])
    return (flo, fhi)

def get_units(ds, vn):
    try:
        units = ds[vn].units
    except AttributeError:
        units = ''
    return units

def add_bathy_contours(ax, ds, depth_levs = [50, 200, 2000], txt=False):
    # This should work with ds being a history file Dataset, or the G dict.
    h = ds.h.values
    # trim depth_levs to be in range of h
    hmax = h.max()
    depth_levs = [item for item in depth_levs if item < hmax]
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    if len(depth_levs) > 0:
        cs = ax.contour(lon, lat, h, depth_levs, colors='k',
            linewidths=0.5) # could add linestyles='dashed' e.g.
        if txt==True:
            ii = 0
            for lev in depth_levs:
                ax.text(.95, .95 - ii*.03, str(lev)+' m', ha='right', transform=ax.transAxes)
                ii += 1
        
def add_map_field(ax, ds, vn, vlims_dict, slev=-1, cmap='rainbow', fac=1,
    alpha=1, aa=[], vlims_fac=3, do_mask_edges=False):
    cmap = plt.get_cmap(name=cmap)
    if 'lon_rho' in ds[vn].coords:
        tag = 'rho'
    if 'lon_u' in ds[vn].coords:
        tag = 'u'
    if 'lon_v' in ds[vn].coords:
        tag = 'v'
        
    x = ds['lon_'+tag].values
    y = ds['lat_'+tag].values
    px, py = get_plon_plat(x,y)
    m = ds['mask_'+tag].values
    if vn in ['zeta', 'ubar', 'vbar']:
        v = ds[vn][0,:,:].values
    else:
        v = ds[vn][0, slev,:,:].values
    v_scaled = fac*v
    
    # account for WET_DRY
    if (tag=='rho') and ('wetdry_mask_rho' in ds.data_vars):
        mwd = ds.wetdry_mask_rho[0,:,:].values.squeeze()
        v_scaled[mwd==0] = np.nan
    
    # SETTING COLOR LIMITS
    # First see if they are already set. If so then we are done.
    vlims = vlims_dict[vn]
    if len(vlims) == 0:
        # If they are not set then set them.
        if len(aa) == 4:
            # make a mask to isolate field for chosing color limits
            x0 = aa[0]; x1 = aa[1]
            y0 = aa[2]; y1 = aa[3]
            m[x<x0] = 0; m[x>x1] = 0
            m[y<y0] = 0; m[y>y1] = 0
            # set section color limits
            fldm = v_scaled[m[1:,1:]==1]
            vlims = auto_lims(fldm, vlims_fac=vlims_fac)
        else:
            vlims = auto_lims(v_scaled, vlims_fac=vlims_fac)
        vlims_dict[vn] = vlims
        # dicts have essentially global scope, so setting it here sets it everywhere
                
    if do_mask_edges:
        v_scaled = mask_edges(v_scaled, x, y)
    
        
    
    cs = ax.pcolormesh(px, py, v_scaled, vmin=vlims[0], vmax=vlims[1], cmap=cmap, alpha=alpha)
    return cs

def add_velocity_vectors(ax, ds, fn, v_scl=10, v_leglen=0.5, nngrid=80, zlev='top', center=(.8,.05)):
    # v_scl: scale velocity vector (smaller to get longer arrows)
    # v_leglen: m/s for velocity vector legend
    xc = center[0]
    yc = center[1]
    # GET DATA
    G = zrfun.get_basic_info(fn, only_G=True)
    if zlev == 'top':
        u = ds['u'][0, -1, :, :].values
        v = ds['v'][0, -1, :, :].values
    elif zlev == 'bot':
        u = ds['u'][0, 0, :, :].values
        v = ds['v'][0, 0, :, :].values
    else:
        zfull_u = get_zfull(ds, fn, 'u')
        zfull_v = get_zfull(ds, fn, 'v')
        u = get_laym(ds, zfull_u, ds['mask_u'].values, 'u', zlev)
        v = get_laym(ds, zfull_v, ds['mask_v'].values, 'v', zlev)
    # ADD VELOCITY VECTORS
    # set masked values to 0
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    # create interpolant
    import scipy.interpolate as intp
    ui = intp.interp2d(G['lon_u'][0, :], G['lat_u'][:, 0], u)
    vi = intp.interp2d(G['lon_v'][0, :], G['lat_v'][:, 0], v)
    # create regular grid
    aaa = ax.axis()
    daax = aaa[1] - aaa[0]
    daay = aaa[3] - aaa[2]
    axrat = np.cos(np.deg2rad(aaa[2])) * daax / daay
    x = np.linspace(aaa[0], aaa[1], int(round(nngrid * axrat)))
    y = np.linspace(aaa[2], aaa[3], int(nngrid))
    xx, yy = np.meshgrid(x, y)
    # interpolate to regular grid
    uu = ui(x, y)
    vv = vi(x, y)
    mask = uu != 0
    # plot velocity vectors
    Q = ax.quiver(xx[mask], yy[mask], uu[mask], vv[mask],
        units='width', scale=v_scl, scale_units='width', color='m')
    plt.quiverkey(Q, .9, .1, v_leglen, str(v_leglen)+' $ms^{-1}$', angle=20)

def add_windstress_flower(ax, ds, t_scl=1, t_leglen=0.1, center=(.85,.25), fs=12):
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
    tauxm = ds.sustr.mean()
    tauym = ds.svstr.mean()
    x = center[0]
    y = center[1]
    ax.quiver([x, x] , [y, y], [tauxm, tauxm], [tauym, tauym],
        units='width', scale=t_scl, scale_units='width', color='k',
        transform=ax.transAxes)
    tt = 1./np.sqrt(2)
    t_alpha = 0.4
    ax.quiver(x*np.ones(8) , y*np.ones(8),
        t_leglen*np.array([0,tt,1,tt,0,-tt,-1,-tt]),
        t_leglen*np.array([1,tt,0,-tt,-1,-tt,0,tt]),
        units='width', scale=t_scl, scale_units='width', color='k', alpha=t_alpha,
        transform=ax.transAxes)
    ax.text(x, y-.13,'Windstress',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes, fontsize=fs)
    ax.text(x, y-.1, str(t_leglen) + ' Pa',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes, fontsize=fs)

def add_info(ax, fn, fs=12, loc='lower_right', his_num=False):
    # put info on plot
    T = zrfun.get_basic_info(fn, only_T=True)
    dt_local = get_dt_local(T['dt'])
    if loc == 'lower_right':
        ax.text(.95, .075, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.text(.95, .065, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        if his_num:
            ax.text(.95, .01, fn.name.split('.')[0].split('_')[-1],
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes, fontsize=fs,
                bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
            
    elif loc == 'upper_right':
        ax.text(.95, .935, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.text(.95, .925, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(.06, .04, str(fn).split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs,
        bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

def get_dt_local(dt, tzl='US/Pacific'):
    # take a model datetime (assumed to be UTC) and return local datetime
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

def get_aa(ds):
    x = ds['lon_psi'][0,:].values
    y = ds['lat_psi'][:,0].values
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa
    
def get_aa_ex(ds):
    x = ds['lon_psi_ex'][0,:].values
    y = ds['lat_psi_ex'][:,0].values
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa

def get_zfull(ds, fn, which_grid):
    # get zfull field on "which_grid" ('rho', 'u', or 'v')
    G, S, T = zrfun.get_basic_info(fn)
    zeta = 0 * G['h']
    zr_mid = zrfun.get_z(G['h'], zeta, S, only_rho=True)
    zr_bot = -G['h'].reshape(1, G['M'], G['L']).copy()
    zr_top = zeta.reshape(1, G['M'], G['L']).copy()
    zfull0 = make_full((zr_bot, zr_mid, zr_top))
    if which_grid == 'rho':
        zfull = zfull0
    elif which_grid == 'u':
        zfull = zfull0[:, :, 0:-1] + np.diff(zfull0, axis=2)/2
    elif which_grid == 'v':
        zfull = zfull0[:, 0:-1, :] + np.diff(zfull0, axis=1)/2
    return zfull

def get_laym(ds, zfull, mask, vn, zlev):
    # make the layer
    fld_mid = ds[vn].values.squeeze()
    fld = make_full((fld_mid,))
    zlev_a = zlev * np.ones(1)
    lay = get_layer(fld, zfull, zlev_a)
    lay[mask == False] = np.nan
    # Note: if mask came directly from the mask_rho field of a history file it is
    # 1 = water, and 0 = land, so it would be better to add nan's using this as a
    # test, but using mask == False works.
    return lay

def get_layer(fld, zfull, which_z):
    """
    Creates a horizontal slice through a 3D ROMS data field.  It is very fast
    because of the use of "choose"
    Input:
        fld (3D ndarray) of the data field to slice
        z (3D ndarray) of z values (like from make_full)
        which_z (ndarray of length 1) of the z value for the layer
    Output:
        lay (2D ndarray) fld on z == which_z,
            with np.nan where it is not defined
    """
    N, M, L = fld.shape # updates N for full fields
    Nmax = 30
    ii = np.arange(0,N,Nmax)
    ii = np.concatenate((ii,np.array([N,])))
    fld0 = np.nan * np.zeros((M, L), dtype=int)
    fld1 = np.nan * np.zeros((M, L), dtype=int)
    z0 = np.nan * np.zeros((M, L), dtype=int)
    z1 = np.nan * np.zeros((M, L), dtype=int)
    # NOTE: need fewer than 32 layers to use "choose"
    # so we split the operation into steps in this loop
    j = 0
    while j < len(ii)-1:
        i_lo = ii[j]
        i_hi = min(ii[j+1] + 1, ii[-1]) # overlap by 1
        NN = i_hi - i_lo # the number of levels in this chunk
        this_zr = zfull[i_lo:i_hi].copy()
        this_fld = fld[i_lo:i_hi].copy()
        zm = this_zr < which_z
        ind0 = np.zeros((M, L), dtype=int)
        ind1 = np.zeros((M, L), dtype=int)
        ind0 = (zm == True).sum(0) - 1 # index of points below which_z
        ind1 = ind0 + 1 # index of points above which_z
        # dealing with out-of-bounds issues
        # note 0 <= ind0 <= NN-1
        # and  1 <= ind1 <= NN
        # make ind1 = ind0 for out of bounds cases
        ind0[ind0 == -1] = 0 # fix bottom case
        ind1[ind1 == NN] = NN-1 # fix top case
        # and now cells that should be masked have equal indices
        this_mask = ind0 != ind1
        this_fld0 = ind0.choose(this_fld)
        this_fld1 = ind1.choose(this_fld)
        this_z0 = ind0.choose(this_zr)
        this_z1 = ind1.choose(this_zr)
        fld0[this_mask] = this_fld0[this_mask]
        fld1[this_mask] = this_fld1[this_mask]
        z0[this_mask] = this_z0[this_mask]
        z1[this_mask] = this_z1[this_mask]
        j += 1
    # do the interpolation
    dz = z1 - z0
    dzf = which_z - z0
    dz[dz == 0] = np.nan
    fr = dzf / dz
    lay = fld0*(1 - fr) + fld1*fr
    return lay

def make_full(flt):
    """
    Adds top and bottom layers to array fld. This is intended for 3D ROMS data
    fields that are on the vertical rho grid, and where we want (typically for
    plotting purposes) to extend this in a smart way to the sea floor and the
    sea surface.
    NOTE: input is always a tuple.  If just sending a single array pack it
    as zfun.make_full((arr,))
    Input:
        flt is a tuple with either 1 ndarray (fld_mid,),
        or 3 ndarrays (fld_bot, fld_mid, fld_top)
    Output: fld is the "full" field
    """
    if len(flt)==3:
       fld = np.concatenate(flt, axis=0)
    elif len(flt)==1:
        if len(flt[0].shape) == 3:
            fld_mid = flt[0]
            N, M, L = fld_mid.shape
            fld_bot = fld_mid[0].copy()
            fld_bot = fld_bot.reshape(1, M, L).copy()
            fld_top = fld_mid[-1].copy()
            fld_top = fld_top.reshape(1, M, L).copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
        elif len(flt[0].shape) == 2:
            fld_mid = flt[0]
            N, M = fld_mid.shape
            fld_bot = fld_mid[0].copy()
            fld_bot = fld_bot.reshape(1, M).copy()
            fld_top = fld_mid[-1].copy()
            fld_top = fld_top.reshape(1, M).copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
        elif len(flt[0].shape) == 1:
            fld_mid = flt[0]
            fld_bot = fld_mid[0].copy()
            fld_top = fld_mid[-1].copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
    return fld
    
def mask_edges(fld, lon, lat):
    """
    Mask out map fields at open ocean edges.
    WARNING: currently hard-coded to only mask WSN edges.
    Input:
        2D fields of data (masked array), and associated lon and lat
        all must be the same shape
    Output:
        The data field, now masked at edges.
    """
    mm = 6
    mask = (lon < lon[0,mm]) | (lat < lat[mm,0]) | (lat > lat[-mm,0])
    fld[mask] = np.nan
    return fld
    
def get_sect(fn, vn, x_e, y_e):
    """
    Code to extract a section of any variable from a history file and return
    fields useful for plotting.
    
    The section plotting concept is that we form fields:
        dist_se: array of distance along section [km] (layer, dist) on cell edges
        zw_se: array of z position along section [m] (layer, dist) on cell edges
        fld_s: array of field values in section (layer, dist), on cell centers
    These are then what we use in pcolormesh.

    Naming conventions:
        [] is a vector on cell centers
        []_e is a vector on cell edges
        []_s is an array on cell centers
        []_se is an array on cell edges

    This should work for variables on any of the horizontal grids (u, v, rho).
    """
    import xarray as xr
    G, S, T = zrfun.get_basic_info(fn)
    ds = xr.open_dataset(fn)

    # Make points x_e, y_e between for interpolating the field onto, on cell centers.
    x = x_e[:-1] + np.diff(x_e)/2
    y = y_e[:-1] + np.diff(y_e)/2

    # Gather some fields, making sure we use the appropriate lon, lat grids.
    if 'eta_u' in ds[vn].dims:
        lon = G['lon_u']
        lat = G['lat_u']
    elif 'eta_v' in ds[vn].dims:
        lon = G['lon_v']
        lat = G['lat_v']
    elif 'eta_rho' in ds[vn].dims:
        lon = G['lon_rho']
        lat = G['lat_rho']
    mask = G['mask_rho']
    h = G['h']
    zeta = ds['zeta'].values.squeeze()
    # Do this to make sure we don't end up with nan's in our zw field.
    h[mask==0] = 0
    zeta[mask==0] = 0

    # get zw field, used for cell edges.
    zw = zrfun.get_z(h, zeta, S, only_w=True)
    N = zw.shape[0]

    # Get 3-D field that we will full the section from
    fld = ds[vn].values.squeeze()
    # Force fld to be on the s_rho grid in the vertical if it is not already.
    if 's_w' in ds[vn].dims:
        fld = fld[:-1,:,:] + np.diff(fld,axis=0)/2

    def get_dist(x,y):
        # Create a vector of distance [km] along a track
        # defined by lon, lat points (x and y in the arguments)
        earth_rad = zfun.earth_rad(np.mean(y)) # m
        xrad = np.pi * x /180
        yrad = np.pi * y / 180
        dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
        dy = earth_rad * np.diff(yrad)
        ddist = np.sqrt(dx**2 + dy**2)
        dist = np.zeros(len(x))
        dist[1:] = ddist.cumsum()/1000 # km
        return dist
    
    dist_e = get_dist(x_e,y_e) # cell edges
    dist = get_dist(x,y) # cell centers

    # Make dist_e into a 2-D array for plotting.
    dist_se = np.ones((N,1)) * dist_e.reshape((1,-1))
    # the -1 means infer size from the array

    def get_sect(x, y, fld, lon, lat):
        # Interpolate a 2-D or 3-D field along a 2-D track
        # defined by lon and lat vectors x and y.
        # We assume that the lon, lat arrays are plaid.
        col0, col1, colf = zfun.get_interpolant(x, lon[1,:])
        row0, row1, rowf = zfun.get_interpolant(y, lat[:,1])
        colff = 1 - colf
        rowff = 1 - rowf
        if len(fld.shape) == 3:
            fld_s = (rowff*(colff*fld[:, row0, col0] + colf*fld[:, row0, col1])
                + rowf*(colff*fld[:, row1, col0] + colf*fld[:, row1, col1]))
        elif len(fld.shape) == 2:
            fld_s = (rowff*(colff*fld[row0, col0] + colf*fld[row0, col1])
                + rowf*(colff*fld[row1, col0] + colf*fld[row1, col1]))
        return fld_s

    # Do the section extractions for zw (edges) and sv (centers)
    zw_se = get_sect(x_e, y_e, zw, lon, lat)
    fld_s = get_sect(x, y, fld, lon, lat )

    # Also generate top and bottom lines, with appropriate masking
    zbot = -get_sect(x, y, h, lon, lat )
    ztop = get_sect(x, y, zeta, lon, lat )
    zbot[np.isnan(fld_s[-1,:])] = 0
    ztop[np.isnan(fld_s[-1,:])] = 0
    
    return x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat
    
def get_section(ds, vn, x, y, in_dict):
    """
    This is OBSOLETE, and has been replaced by the much cleaner get_sect() above.
    """

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    h = G['h']
    zeta = ds['zeta'].values.squeeze()
    zr = zrfun.get_z(h, zeta, S, only_rho=True)

    sectvar = ds[vn].values.squeeze()

    L = G['L']
    M = G['M']
    N = S['N']

    lon = G['lon_rho']
    lat = G['lat_rho']
    mask = G['mask_rho']
    maskr = mask.reshape(1, M, L).copy()
    mask3 = np.tile(maskr, [N, 1, 1])
    zbot = -h # don't need .copy() because of the minus operation

    # make sure fields are masked
    zeta[mask==0] = np.nan
    zbot[mask==0] = np.nan
    sectvar[mask3==0] = np.nan

    # create dist
    earth_rad = zfun.earth_rad(np.mean(lat[:,0])) # m
    xrad = np.pi * x /180
    yrad = np.pi * y / 180
    dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
    dy = earth_rad * np.diff(yrad)
    ddist = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(len(x))
    dist[1:] = ddist.cumsum()/1000 # km
    # find the index of zero
    i0, i1, fr = zfun.get_interpolant(np.zeros(1), dist)
    idist0 = i0
    distr = dist.reshape(1, len(dist)).copy()
    dista = np.tile(distr, [N, 1]) # array
    # pack fields to process in dicts
    d2 = dict()
    d2['zbot'] = zbot
    d2['zeta'] = zeta
    d2['lon'] = lon
    d2['lat'] = lat
    d3 = dict()
    d3['zr'] = zr
    d3['sectvar'] = sectvar
    # get vectors describing the (plaid) grid
    xx = lon[1,:]
    yy = lat[:,1]
    col0, col1, colf = zfun.get_interpolant(x, xx)
    row0, row1, rowf = zfun.get_interpolant(y, yy)
    # and prepare them to do the bilinear interpolation
    colff = 1 - colf
    rowff = 1 - rowf
    # now actually do the interpolation
    # 2-D fields
    v2 = dict()
    for fname in d2.keys():
        fld = d2[fname]
        fldi = (rowff*(colff*fld[row0, col0] + colf*fld[row0, col1])
        + rowf*(colff*fld[row1, col0] + colf*fld[row1, col1]))
        v2[fname] = fldi
    # 3-D fields
    v3 = dict()
    for fname in d3.keys():
        fld = d3[fname]
        fldid = (rowff*(colff*fld[:, row0, col0] + colf*fld[:, row0, col1])
        + rowf*(colff*fld[:, row1, col0] + colf*fld[:, row1, col1]))
        v3[fname] = fldid
    v3['dist'] = dista # distance in km
    # make "full" fields by padding top and bottom
    nana = np.nan * np.ones((N + 2, len(dist))) # blank array
    v3['zrf'] = nana.copy()
    v3['zrf'][0,:] = v2['zbot']
    v3['zrf'][1:-1,:] = v3['zr']
    v3['zrf'][-1,:] = v2['zeta']
    #
    v3['sectvarf'] = nana.copy()
    v3['sectvarf'][0,:] = v3['sectvar'][0,:]
    v3['sectvarf'][1:-1,:] = v3['sectvar']
    v3['sectvarf'][-1,:] = v3['sectvar'][-1,:]
    #
    v3['distf'] = nana.copy()
    v3['distf'][0,:] = v3['dist'][0,:]
    v3['distf'][1:-1,:] = v3['dist']
    v3['distf'][-1,:] = v3['dist'][-1,:]    
    
    # attempt to skip over nan's
    v3.pop('zr')
    v3.pop('sectvar')
    v3.pop('dist')
    mask3 = ~np.isnan(v3['sectvarf'][:])
    #print(mask3.shape)
    mask2 = mask3[-1,:]
    dist = dist[mask2]
    NC = len(dist)
    NR = mask3.shape[0]
    for k in v2.keys():
        #print('v2 key: ' + k)
        v2[k] = v2[k][mask2]
    for k in v3.keys():
        #print('v3 key: ' + k)
        v3[k] = v3[k][mask3]
        v3[k] = v3[k].reshape((NR, NC))
        #print(v3[k].shape)
    
    return v2, v3, dist, idist0
    
def maxmin(a):
    # find the value and location of the max and min of a 2D
    # masked array
    try:
        jmax,imax = np.unravel_index(np.nanargmax(a),a.shape)
        amax = a[jmax,imax]
        jmin,imin = np.unravel_index(np.nanargmin(a),a.shape)
        amin = a[jmin,imin]
    except ValueError:
        amax, jmax, imax, amin, jmin, imin = (0,0,0,0,0,0)
    return amax, jmax, imax, amin, jmin, imin

def draw_box(ax, aa, linestyle='-', color='k', alpha=1, linewidth=.5, inset=0):
    aa = [aa[0]+inset, aa[1]-inset, aa[2]+inset, aa[3]-inset]
    ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]], [aa[2], aa[2], aa[3], aa[3], aa[2]],
        linestyle=linestyle, color=color, alpha=alpha, linewidth=linewidth)
        


