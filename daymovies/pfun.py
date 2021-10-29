"""
Module of basic utilities for plotting.  The goal is to make the code in
pfun less cumbersome to write and edit.

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun
import netCDF4 as nc

import numpy as np
from datetime import datetime, timedelta
import pytz

import ephem_functions as efun

if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # regular (remote, linux) version
    import matplotlib as mpl
    mpl.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import pinfo

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

def get_tracks(Q, Ldir):
    if Q['test'] == False:
        dt0 = datetime.strptime(Q['ds0'],'%Y.%m.%d')
        dt1 = datetime.strptime(Q['ds1'],'%Y.%m.%d')
        dtt = (dt1-dt0).days + 1
        import subprocess
        cmd = ['python', '../tracker/tracker.py', '-exp', Q['exp'],
            '-ds', Q['ds0'], '-dtt', str(dtt), '-t', Q['ttag'], '-clb', 'True']
        #print(cmd)
        proc = subprocess.Popen(cmd)
        proc.communicate()
    else:
        # using test = True allows us to use pre-calculated tracks
        pass
    tr_fn = Ldir['LOo'] + 'tracks/' + Q['exp'] + '_surf_' + Q['ttag'] + '/release_' + Q['ds0'] + '.nc'
    #print(tr_fn)
    Q['tr_fn'] = tr_fn
    
def get_speed(ds, nlev):
    u = ds['u'][0,nlev,:,:]
    v = ds['v'][0,nlev,:,:]
    u[u.mask] = 0
    v[v.mask] = 0
    fld0 = 0 * ds['salt'][0,-1,1:-1,1:-1]
    uu = (u[1:-1,:-1] + u[1:-1,1:])/2
    vv = (v[:-1,1:-1] + v[1:,1:-1])/2
    fld = np.sqrt(uu*uu + vv*vv)
    fld = np.ma.masked_where(fld0.mask, fld)
    return fld

def get_arag(ds, Q, aa, nlev):
    from PyCO2SYS import CO2SYS
    import seawater as sw
    G = zrfun.get_basic_info(Q['fn'], only_G=True)
    # find indices that encompass region aa
    i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
    i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
    j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
    j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
    px = G['lon_psi'][j0:j1-1, i0:i1-1]
    py = G['lat_psi'][j0:j1-1, i0:i1-1]
    lat = G['lat_rho'][j0:j1,i0:i1] # used in sw.pres
    # first extract needed fields and save in v_dict
    v_dict = {}
    vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
    for cvn in vn_in_list:
        L = zfun.fillit(ds[cvn][0,nlev,j0:j1,i0:i1])
        v_dict[cvn] = L
    # ------------- the CO2SYS steps -------------------------
    # create pressure
    Ld = G['h'][j0:j1,i0:i1]
    Lpres = sw.pres(Ld, lat)
    # get in situ temperature from potential temperature
    Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
    # convert from umol/L to umol/kg using in situ dentity
    Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
    Lalkalinity[Lalkalinity < 100] = np.nan
    LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
    LTIC[LTIC < 100] = np.nan
    Lpres = zfun.fillit(Lpres)
    Ltemp = zfun.fillit(Ltemp)
    CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
        Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
    # PH = CO2dict['pHout']
    # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
    ARAG = CO2dict['OmegaARout']
    ARAG = zfun.fillit(ARAG.reshape((v_dict['salt'].shape)))
    fld = ARAG[1:-1, 1:-1]
    return px, py, fld

def get_ax_limits(Q):
    # Set limits and ticklabels, and other info.
    # A good aspect ratio is dlat/dlon = 1/.8
    if Q['dom'] == 'full':
        Q['aa'] = []
        Q['xtl'] = range(-129,-121,2)
        Q['ytl'] = range(44,52,2)
        Q['v_scl'] = 3 # used for velocity vectors
        Q['exp'] = 'p5_merhab' # particle tracking
    elif Q['dom'] == 'PS':
        Q['aa'] = [-123.6, -122, 47, 49]
        Q['xtl'] = [-123, -122.5]
        Q['ytl'] = [47.5, 48, 48.5]
        Q['v_scl'] = 25
        Q['exp'] = 'p5_PS'
    elif Q['dom'] == 'Psouth':
        # South Sound
        Q['aa'] = [-123.15, -122.5, 47, 47.5]#.8125]
        Q['xtl'] = [-123, -122.5]
        Q['ytl'] = [47, 47.5]
        Q['v_scl'] = 40
    elif Q['dom'] == 'willapa':
        Q['aa'] = [-124.6, -123.65, 46, 47.2]
        Q['xtl'] = [-124.4, -124.2, -124, -123.8]
        Q['ytl'] = [46, 46.5, 47]
        Q['v_scl'] = 25
    elif Q['dom'] == 'nshelf':
        # north part of the Washington shelf
        Q['aa'] = [-126.1, -123.7, 45.8, 48.8]
        Q['xtl'] = [-126, -125, -124]
        Q['ytl'] = [46, 47, 48]
        Q['v_scl'] = 25
        
def get_moor_info(Q):
    # set mooring info
    M = dict()
    if Q['dom'] == 'full':
        M['lon'] = -124.5
        M['lat'] = 47
        M['city'] = 'Westport'
        M['wscl'] = 20
    elif Q['dom'] =='PS':
        M['lon'] = -122.433
        M['lat'] = 47.86
        M['city'] = 'Seattle'
        M['wscl'] = 10
    elif Q['dom'] =='Psouth':
        M['lon'] = -122.712
        M['lat'] = 47.31
        M['city'] = 'Seattle'
        M['wscl'] = 10
    elif Q['dom'] in ['willapa', 'nshelf']:
        M['lon'] = -124.4
        M['lat'] = 46.6
        M['city'] = 'Westport'
        M['wscl'] = 20
    return M

def plot_time_series(ax, M, T):
    iot = zfun.find_nearest_ind(M['ot'], T['ocean_time'])
    zeta = M['zeta']*3.28084 # convert meters to feet
    ax.plot(M['ot'], zeta, '-k', lw=2)
    ax.plot(T['ocean_time'], zeta[iot],'ok', ms=20)
    x0 = M['ot'][0]; x1 = M['ot'][-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta.min())-5
    y1 = np.ceil(zeta.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], 'k', alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(M['ot'][0] + (M['ot'][-1]-M['ot'][0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(M['ot'][0] + (M['ot'][-1]-M['ot'][0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    dt_local = get_dt_local(T['tm'])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p')+' '+dt_local.tzname(),
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c='k',
        transform=ax.transAxes)
    # add day/night
    Srise = M['Srise']
    Sset = M['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

def get_moor(ds0, ds1, Ldir, Q, M):
    """
    Gets a simple mooring record somewhere in the domain, to use in the plotting.
    It packs all its data in the dictionary Q, which - because of the odd
    scope behaviour of python dictionaries, is modified even as far as the
    calling function is concerned, and so does not have to be returned.
    """
    m_fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
    ot_list = []
    dt_list = []
    zeta_list = []
    uwind_list = []
    vwind_list = []
    G = zrfun.get_basic_info(m_fn_list[0], only_G=True)
    
    mi = zfun.find_nearest_ind(G['lon_rho'][0,:], M['lon'])
    mj = zfun.find_nearest_ind(G['lat_rho'][:,0], M['lat'])
    for fn in m_fn_list:
        T = zrfun.get_basic_info(fn, only_T=True)
        dt_list.append(T['tm'])
        ds = nc.Dataset(fn)
        ot_list.append(ds['ocean_time'][0])
        zeta_list.append(ds['zeta'][0,mj,mi])
        uwind_list.append(ds['Uwind'][0,mj,mi])
        vwind_list.append(ds['Vwind'][0,mj,mi])
        ds.close()
    ot = zfun.fillit(np.array(ot_list))
    zeta = zfun.fillit(np.array(zeta_list))
    uwind = zfun.fillit(np.array(uwind_list))
    vwind = zfun.fillit(np.array(vwind_list))
    uwind = zfun.filt_hanning(uwind, n=5, nanpad=False)
    vwind = zfun.filt_hanning(vwind, n=5, nanpad=False)
    M['ot'] = ot
    M['zeta'] = zeta
    M['uwind'] = uwind
    M['vwind'] = vwind
    
    # get sunrise and sunset
    sdt0 = dt_list[0] - timedelta(days=1)
    sdt1 = dt_list[-1] + timedelta(days=1)
    Srise, Sset = efun.get_sunrise_sunset(sdt0, sdt1, city=M['city'])
    # these are lists of datetimes in the UTC time zone
    M['Srise'] = Srise
    M['Sset'] = Sset

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.cos(np.pi*yav/180))

def add_coast(ax, dir0=Ldir['data'], color='k', lw=0.5):
    fn = dir0 + 'coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    ax.plot(C['lon'].values, C['lat'].values, '-', color=color, linewidth=lw)

def mask_edges(ds, fld, Q):
    # mask off selected edges, e.g. for NPZD variables
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    aa = Q['aa']
    clat = np.cos((np.pi/180)*(aa[2]+aa[3])/2)
    pad = .3
    fld = np.ma.masked_where(xr<aa[0]+pad, fld)
    #fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2]+pad*clat, fld)
    fld = np.ma.masked_where(yr>aa[3]-pad*clat, fld)
    return fld

def get_vlims(ds, fld, Q):
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    aa = Q['aa']
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    (vmin, vmax) = auto_lims(fld, vlims_fac=pinfo.vlims_fac_dict[Q['vn']])
    Q['vmin'] = vmin
    Q['vmax'] = vmax

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

def add_bathy_contours(ax, ds, depth_levs = [200], txt=False, c='k', lw=0.5):
    # this should work with ds being a history file Dataset, or the G dict.
    h = ds['h'][:]
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    c = c
    ax.contour(lon, lat, h, depth_levs, colors=c, linewidths=lw)
    if txt==True:
        ii = 0
        for lev in depth_levs:
            ax.text(.95, .95 - ii*.03, str(lev)+' m', c=c,
                    ha='right', transform=ax.transAxes)
            ii += 1

def add_velocity_vectors(ax, aa, ds, fn, v_scl=3, nngrid=80, zlev='top'):
    # v_scl: scale velocity vector (smaller to get longer arrows)
    # GET DATA
    G = zrfun.get_basic_info(fn, only_G=True)
    if zlev == 'top':
        u = ds['u'][0, -1, :, :].squeeze()
        v = ds['v'][0, -1, :, :].squeeze()
    elif zlev == 'bot':
        u = ds['u'][0, 0, :, :].squeeze()
        v = ds['v'][0, 0, :, :].squeeze()
    # ADD VELOCITY VECTORS
    # set masked values to 0
    ud = u.data; ud[u.mask]=0
    vd = v.data; vd[v.mask]=0
    # create interpolant
    import scipy.interpolate as intp
    ui = intp.interp2d(G['lon_u'][0, :], G['lat_u'][:, 0], ud)
    vi = intp.interp2d(G['lon_v'][0, :], G['lat_v'][:, 0], vd)
    # create regular grid
    daax = aa[1] - aa[0]
    daay = aa[3] - aa[2]
    axrat = np.cos(np.deg2rad(aa[2])) * daax / daay
    x = np.linspace(aa[0], aa[1], int(round(nngrid * axrat)))
    y = np.linspace(aa[2], aa[3], int(nngrid))
    xx, yy = np.meshgrid(x, y)
    # interpolate to regular grid
    uu = ui(x, y)
    vv = vi(x, y)
    mask = uu != 0
    # plot velocity vectors
    ax.quiver(xx[mask], yy[mask], uu[mask], vv[mask],
        units='y', scale=v_scl, scale_units='y', color='b')
        
def add_wind(ax, M, T):
    """Add a windspeed vector with circles for scale."""
    # scl is windspeed [knots] for a 1 inch circle or arrow
    # this makes a circle 1 inch (72 points) in radius
    # and a vector 1 inch long for a windspeed of "scl" knots
    iot = zfun.find_nearest_ind(M['ot'], T['ocean_time'])
    uwind = M['uwind'][iot]
    vwind = M['vwind'][iot]
    ax.plot(M['lon'],M['lat'],'o', ms=144, mfc='None', mec='k', mew=1.5, alpha=.6)
    ax.quiver(M['lon'],M['lat'], uwind*1.94384, vwind*1.94384,
            scale=M['wscl'], scale_units='inches',
            headwidth=5,headlength=5, color='k')
            
def add_wind_text(ax, aa, M, fs):
    x0 = M['lon']
    y0 = M['lat']
    dx = aa[3] - aa[2]
    xt = x0 + dx/9
    yt = y0 - dx/30
    ax.text(xt, yt, str(M['wscl'])+' knot\nwind', c='k',
        ha='center', va='center', style='italic', size=.7*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=.3))

def add_info(ax, fn, fs=12, loc='lower_right'):
    # put info on plot
    T = zrfun.get_basic_info(fn, only_T=True)
    dt_local = get_dt_local(T['tm'])
    if loc == 'lower_right':
        ax.text(.95, .075, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .065, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
    elif loc == 'upper_right':
        ax.text(.95, .935, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .925, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)

def get_dt_local(dt, tzl='US/Pacific'):
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

def get_aa(ds):
    x = ds['lon_psi'][0,:]
    y = ds['lat_psi'][:,0]
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa

def draw_box(ax, aa, linestyle='-', color='k', alpha=1, linewidth=.5, inset=0):
    aa = [aa[0]+inset, aa[1]-inset, aa[2]+inset, aa[3]-inset]
    ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]], [aa[2], aa[2], aa[3], aa[3], aa[2]],
        linestyle=linestyle, color=color, alpha=alpha, linewidth=linewidth)