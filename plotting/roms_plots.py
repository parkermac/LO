"""
Module of plotting functions.

Each function creates, and optionally saves, a plot of fields
from a ROMS history file.

INPUT: in_dict: a tuple with information to pass to the plot, such as:
- fn: text string with the full path name of the history file to plot
- fn_out: text string with full path of output file name
- auto_vlims: a boolean governing how color limits are set
- testing: a boolean for testing (e.g. shorter, faster particle tracking)
OUTPUT: either a screen image or a graphics file

"""
import numpy as np
import xarray as xr
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
import pinfo
from importlib import reload
reload(pfun)
reload(pinfo)

Ldir = Lfun.Lstart()
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
    
def P_basic(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['salt', 'temp']
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_ri(in_dict):
    """
    Simplified Richardson number
    """
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(20,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PLOT CODE
    xrho = ds['lon_rho'][0,:].values
    yrho = ds['lat_rho'][:,0].values

    # define box
    aa = [-123.25, -122.1, 47, 48.75]
    ix0 = zfun.find_nearest_ind(xrho, aa[0])
    ix1 = zfun.find_nearest_ind(xrho, aa[1])
    iy0 = zfun.find_nearest_ind(yrho, aa[2])
    iy1 = zfun.find_nearest_ind(yrho, aa[3])

    h = ds.h[iy0:iy1, ix0:ix1].values
    rho_bot = ds.rho[0, 0, iy0:iy1, ix0:ix1].values
    rho_top = ds.rho[0, -1, iy0:iy1, ix0:ix1].values
    drho = rho_bot - rho_top
    u = ds.ubar[0, iy0:iy1, ix0-1:ix1].values
    v = ds.vbar[0, iy0-1:iy1, ix0:ix1].values
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    uu = (u[:, 1:] + u[:, :-1])/2
    vv = (v[1:, :] + v[:-1, :])/2
    spd2 = uu**2 + vv**2
    spd2[np.isnan(drho)] = np.nan
    spd2[spd2 < .001] = .001 # avoid divide by zero errors

    # approximate Richardson number
    rho0 = ds.rho0.values
    g = 9.8
    Ri = g * drho * h / (rho0 * spd2)

    # psi_grid coordinates
    x, y = np.meshgrid(ds.lon_u.values[0,ix0-1:ix1], ds.lat_v.values[iy0-1:iy1,0])

    # PLOTTING
    plt.close('all')
    pfun.start_plot(fs=10, figsize=(18,10))
    fig = plt.figure()

    xt = [-123.2, -122.2]
    yt = [47, 47.5, 48, 48.5]

    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(x, y, drho, vmin=0, vmax=5, cmap=cm.dense)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$\Delta\rho\ [kg\ m^{-3}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(x, y, np.sqrt(spd2), vmin=0, vmax=2, cmap=cm.speed)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'Speed $[m\ s^{-1}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])

    ax = fig.add_subplot(133)
    cs = ax.pcolormesh(x, y, 4*Ri, vmin=0, vmax = 2, cmap='RdYlBu')
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$4 x Ri$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])
        
    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_Chl_DO(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    fs = 14
    ii = 1
    for vn in vn_list:
        if vn == 'phytoplankton':
            slev = -1
            stext = 'Surface'
        elif vn == 'oxygen':
            slev = 0
            stext = 'Bottom'
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        pfun.add_bathy_contours(ax, ds, txt=True)
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_ths(in_dict):
    # Plot property-property plots, like theta vs. s
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    # make a potential density field
    import seawater as sw
    s0 = 25; s1 = 35
    th0 = 0; th1 = 20
    SS, TH = np.meshgrid(np.linspace(s0, s1, 50), np.linspace(th0, th1, 50))
    SIG = sw.dens0(SS, TH) - 1000
    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    h = ds['h'].values
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    s = ds['salt'].values.squeeze()
    th = ds['temp'].values.squeeze()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Theta (deg C)')
    ax.contour(SS, TH, SIG, 20)
    nsub = 500
    alpha = .1
    mask = z > -10
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.r', alpha=alpha)
    mask = (z < -10) & (z > -200)
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.g', alpha=alpha)
    mask = z < -200
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.b', alpha=alpha)
    ax.set_xlim(s0, s1)
    ax.set_ylim(th0, th1)
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_debug(in_dict):
    # Focused on debugging

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=pinfo.figsize)
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
            # add debugging information to the plot
            # gathering info
            u = ds['u'][0,-1,:,:].values
            umax, ujmax, uimax, umin, ujmin, uimin = pfun.maxmin(u)
            v = ds['v'][0,-1,:,:].values
            vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
            eta = ds['zeta'][0,:,:].values
            emax, ejmax, eimax, emin, ejmin, eimin = pfun.maxmin(eta)
            #
            G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
            def add_info(G, ax, name, grd, vval, vj, vi, ypos, clr):
                ax.text(.98, ypos,'%s = %5.1f' % (name, vval), fontweight='bold',
                    horizontalalignment='right', transform=ax.transAxes, color=clr)
                ax.plot(G['lon_'+grd][vj,vi],G['lat_'+grd][vj,vi],'*', color=clr,
                    markeredgecolor='w',markersize=18,)
            add_info(G, ax, 'umax', 'u', umax, ujmax, uimax, .36, 'r')
            add_info(G, ax, 'umin', 'u', umin, ujmin, uimin, .33, 'orange')
            add_info(G, ax, 'vmax', 'v', vmax, vjmax, vimax, .3, 'b')
            add_info(G, ax, 'vmin', 'v', vmin, vjmin, vimin, .27, 'g')
            add_info(G, ax, 'emax', 'rho', emax, ejmax, eimax, .24, 'k')
            add_info(G, ax, 'emin', 'rho', emin, ejmin, eimin, .21, 'm')
        ii += 1

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['oxygen', 'temp']
    z_level = -250
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'phytoplankton'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -3500
        x = np.linspace(lon.min(), lon.max(), 500)
        y = 47 * np.ones(x.shape)
    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] / 'section_lines'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path / track
            # get the track to interpolate onto
            pdict = pickle.load(open(track_fn, 'rb'))
            xx = np.concatenate((xx,pdict['lon_poly']))
            yy = np.concatenate((yy,pdict['lat_poly']))
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    pfun.add_coast(ax)
    aaf = [-125.5, -122.1, 46.8, 50.3] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_splash(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  Eventually I could automate making this new every day.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from PyCO2SYS import CO2SYS
    import seawater as sw
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from warnings import filterwarnings
    filterwarnings('ignore') # skip a warning from PyCO2SYS

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    def get_arag(ds, fn, aa, nlev):
        G = zrfun.get_basic_info(fn, only_G=True)
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
            L = ds[cvn][0,nlev,j0:j1,i0:i1].values
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
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # PH = CO2dict['pHout']
        # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        ARAG = ARAG.reshape((v_dict['salt'].shape))
        ARAG = ARAG[1:-1, 1:-1]
        return px, py, ARAG

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'][0].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.2,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.98,.95,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='right', va='top', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Willapa and Grays Harbor
    aa = [-124.6, -123.65, 46, 47.2]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)

    # SMALL. MAP
    ax = fig.add_subplot(122)
    px, py, ARAG = get_arag(ds, fn, aa, nlev)
    cs = ax.pcolormesh(px,py,ARAG, cmap='coolwarm_r', vmin=0, vmax=3)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-124.5, -124])
    ax.set_yticks([46, 47])
    ax.set_xlabel('Longitude')
    ax.text(.98,.99,'Bottom water\nAragonite\nSaturation\nState',
         ha='right', va='top', weight='bold', transform=ax.transAxes)

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_salt(in_dict):
    # Plot salinity maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'salt'
    vlims = (28.5, 33) # full map
    vlims2 = (22, 32) # PS map
    vlims3 = (29, 32) # PS section
    cmap = 'Spectral_r'

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nSalinity\n'
        + datetime.strftime(T['dt'][0], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'][0], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    
    ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
        va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'][0] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_oxygen(in_dict):
    # Plot bottom oxygen maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'oxygen'
    vlims = (0, 10) # full map
    vlims2 = (0, 10) # PS map
    vlims3 = (0, 10) # PS section
    from cmocean import cm
    cmap = cm.oxy

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_HC_thalweg_long.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, 0, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nBottom Oxygen\n'
        + datetime.strftime(T['dt'][0], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'][0], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=.85*fs)
    ax.text(1, .85, r'$[mg\ L^{-1}]$', transform=ax.transAxes, fontsize=fs, ha='right')
    # ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nHood Canal', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'][0] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_chl(in_dict):
    # Plot phytoplankton maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'phytoplankton'
    vlims = (0, 25) # full map
    vlims2 = (0, 25) # PS map
    vlims3 = (0, 25) # PS section
    cmap = 'Spectral_r'
    
    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nPhytoplankton\n'+pinfo.units_dict[vn]+'\n'
        + datetime.strftime(T['dt'][0], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.99,.97,'range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'][0], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'][0] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
