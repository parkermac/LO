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

# imports

# The calling function, pan_plot.py, has already put alpha on the
# path so it is on the path here too.
import Lfun
Ldir = Lfun.Lstart()
if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # fjord/boiler version
    import matplotlib as mpl
    mpl.use('Agg')
import zfun
import zrfun
from importlib import reload
import pfun; reload(pfun)
import pinfo; reload(pinfo)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm

def P_basic(in_dict):

    # START
    fig = plt.figure(figsize=(18,11)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            vlims_fac=1
        elif vn == 'temp':
            vlims_fac=2.5
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=vlims_fac)
        cb = fig.colorbar(cs)
        cb.ax.tick_params(labelsize=fs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels
        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_AI_Rocky(in_dict):
    """
    This plots maps of surface divergence and vorticity.  It is
    optimized for the old Admiralty Inlet simulations that Dave Sutherland did.
    These don't have the usual LiveOcean naming, but you can plot with commands like:
    
    run pan_plot.py -lt snapshot -pt P_AI_Rocky -g ainlet -t v0 -x old -0 2006.06.01 -hn 700
    
    We use gourad shading, which means that the coordinates have the same size as the data:
    - dive is on the clipped rho grid
    - vort is on the psi grid
    """
    # START
    fs = 16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(14,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa = [-122.8, -122.54, 47.92, 48.22]
    import cmocean
    cmap = cmocean.cm.balance
    # cmap = 'RdYlBu_r'

    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    
    # plot Code
    
    # calculate divergence and vorticity
    uu = ds['u'][0, -1, :, :]
    vv = ds['v'][0, -1, :, :]
    u = zfun.fillit(uu)
    v = zfun.fillit(vv)
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    
    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
    
    dive = ((np.diff(u, axis=1)/G['DX'][:, 1:-1])[1:-1, :]
            + (np.diff(v, axis = 0)/G['DY'][1:-1, :])[:, 1:-1])
    #dive[G['mask_rho'][1:-1,1:-1]==False] = np.nan
    
    vort = np.diff(v, axis=1)/G['DX'][1:,1:] - np.diff(u, axis=0)/G['DY'][1:,1:]
    #vort[G['mask_rho'][1:,1:]==False] = np.nan
    
    scl = 2e-3
    
    # panel 1
    ax = fig.add_subplot(121)
    # cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], dive/scl, cmap=cmap,
    #                     vmin=-1, vmax=1)
    cs = plt.pcolormesh(G['lon_rho'][1:-1,1:-1], G['lat_rho'][1:-1,1:-1], dive/scl, cmap=cmap,
                        vmin=-1, vmax=1, shading='gouraud')
    tstr = (r'Surface Divergence (%0.1e $s^{-1}$)' % (scl))
    #pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])
    ax.set_xticks([-122.8, -122.7, -122.6])
    ax.set_yticks([48, 48.1, 48.2])
    #
    # panel 2
    ax = fig.add_subplot(122)
    # cs = plt.pcolormesh(G['lon_rho'], G['lat_rho'], vort/scl, cmap=cmap,
    #                     vmin=-1, vmax=1)
    cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], vort/scl, cmap=cmap,
                        vmin=-1, vmax=1, shading='gouraud')
    tstr = (r'Surface Vorticity (%0.1e $s^{-1}$)' % (scl))
    ax.set_xticks([-122.8, -122.7, -122.6])
    ax.set_yticks([])
    #fig.colorbar(cs)
    
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    
    #pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr)    
    
    #fig.tight_layout()
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
    
def P_AI_Rocky2(in_dict):
    """
    This plots maps of surface currents and salinity.  It is
    optimized for the old Admiralty Inlet simulations that Dave Sutherland did.
    100 m grid, saves every half hour.  History file 700 is about max ebb at 6:30 PM PDT
    on 2006.06.06.
    
    These don't have the usual LiveOcean naming, but you can plot with commands like:
    
    run pan_plot.py -lt snapshot -pt P_AI_Rocky2 -g ainlet -t v0 -x old -0 2006.06.01 -hn 700
    
    or
    
    run pan_plot.py -lt allhours -pt P_AI_Rocky2 -g ainlet -t v0 -x old -0 2006.06.01 -mov True
    
    """
    # START
    fs = 16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(8,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa = [-122.8, -122.54, 47.92, 48.22]
    # import cmocean
    # cmap = cmocean.cm.speed
    # cmap = 'RdYlBu_r'
    cmap = 'Spectral_r'

    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    
    # plot Code
    
    # calculate speed
    uu = ds['u'][0, -1, :, :]
    vv = ds['v'][0, -1, :, :]
    u = zfun.fillit(uu)
    v = zfun.fillit(vv)
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    # interpolate to the clipped rho grid
    ur = (u[1:-1,1:] + u[1:-1,:-1])/2
    vr = (v[1:,1:-1] + v[:-1,1:-1])/2
    spd = np.sqrt(ur**2 + vr**2)
    spd[spd==0] = np.nan
    
    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
        
    # panel 1
    ax = fig.add_subplot(111)
    cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], spd, cmap=cmap,
                        vmin=0, vmax=1.5)
    tstr = (r'Admiralty Inlet Surface Speed ($m\ s^{-1}$)')
    #pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])
    ax.set_xticks([-122.8, -122.7, -122.6])
    ax.set_yticks([48, 48.1, 48.2])
    
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="4%", height="30%", loc='upper right', borderpad=3)
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    
    #pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr)
    
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=100, v_leglen=1,
        nngrid=80, zlev='top', center=(.1,.05))
    
    fig.tight_layout()
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
    
def P_dye(in_dict):

    # plots a dye field
    # run with -x lo8dye for example

    # START
    fig = plt.figure(figsize=(18,11)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    import cmocean
    cmap = cmocean.cm.haline
    vn_list = ['dye_01', 'dye_01']
    aa = [-123.3, -122.3, 47.2, 48.2] # Hood Canal
    fs = 14
    ii = 1
    for vn in vn_list:
        if ii == 1:
            slev = 0
            ttext = 'Bottom Dye'
        elif ii == 2:
            slev = -1
            ttext = 'Surface Dye'
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=cmap)
        #cb = fig.colorbar(cs)
        #cb.ax.tick_params(labelsize=fs)
        # pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title(ttext, fontsize=fs+2)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels
        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
            # pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            # pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
            pass
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Chl_DO(in_dict):

    # START
    fig = plt.figure(figsize=(18,11)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if ii == 1:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap='ocean_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
            pfun.add_bathy_contours(ax, ds, txt=False)
        elif ii == 2:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    slev=0, cmap='rainbow_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Bottom %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
            pfun.add_bathy_contours(ax, ds, txt=True)
        cb = fig.colorbar(cs)
        cb.ax.tick_params(labelsize=fs)
        pfun.add_coast(ax)
        #ax.axis(pfun.get_aa(ds))
        ax.axis([-129, -122, 42.5, 51.5])
        pfun.dar(ax)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels
        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
            pfun.add_windstress_flower(ax, ds)
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], center=(.8,.85), v_scl=3, v_leglen=1)
        ii += 1
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_DO(in_dict):

    # START
    fs = 18
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(14,12)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa_list = [[-127, -122, 42.5, 51.5], [-123.8, -122.15, 46.85, 49]]
    ii = 1
    import cmocean
    cmap = cmocean.cm.balance_r

    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages

    vn = 'oxygen'
    fac = 0.032
    unit_str = ' $[mg\ L^{-1}]$'
    pinfo.vlims_dict[vn] = (0,4)

    for aa in aa_list:
        vn = 'oxygen'
        ax = fig.add_subplot(1, 2, ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=0, cmap=cmap, fac=fac)
        if ii == 1:
            ax.set_xticks([ -127, -125, -123])
            ax.set_title('(a) Continental Shelf', weight='bold')
            pfun.add_bathy_contours(ax, ds, txt=True)
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
            pfun.add_windstress_flower(ax, ds)
            pfun.draw_box(ax, aa_list[1], color='c', alpha=1, linewidth=3, inset=0)
        elif ii == 2:
            ax.set_title('(b) Puget Sound', color='c', weight='bold')
            ax.set_xticks([-123.5, -123, -122.5])
            ax.set_yticks([47, 48, 49])
            # Inset colorbar
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
            fig.colorbar(cs, cax=cbaxes, orientation='vertical')
            ax.text(.15, .3, 'Bottom %s\n  %s' % (pinfo.tstr_dict[vn],unit_str),
                transform=ax.transAxes, weight='bold')
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ii += 1
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()

def P_rho(in_dict):
    # Surface density and stratification, for Jim Moum 2019_05

    # START
    fig = plt.figure(figsize=(12,8)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa = [-127, -123, 43, 48]
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    rho_surf = ds['rho'][0,-1,:,:] # surface density anomaly
    z_level = -30
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    rho_deep = pfun.get_laym(ds, zfull, ds['mask_rho'][:], 'rho', z_level)
    strat = rho_deep - rho_surf

    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x, y, rho_surf[1:-1,1:-1], vmin=22, vmax=26, cmap='rainbow')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface Density Anomaly (kg m-3)')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds, t_scl=0.3, t_leglen=0.1, center=(0.25, 0.85))

    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x, y, strat[1:-1,1:-1], vmin=0, vmax=3, cmap='jet')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=False)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Density Diff. over top ' + str(-z_level) + ' m (kg m-3)')
    ax.set_xlabel('Longitude')

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_ths(in_dict):
    # Plot property-property plots, like theta vs. s

    # START
    fig = plt.figure(figsize=(12,8)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE

    # make a potential density field
    import seawater as sw
    s0 = 27; s1 = 35
    th0 = 4; th1 = 20
    SS, TH = np.meshgrid(np.linspace(s0, s1, 50), np.linspace(th0, th1, 50))
    SIG = sw.dens0(SS, TH) - 1000

    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    h = ds['h'][:]
    z = zrfun.get_z(h, 0*h, S, only_rho=True)

    s = ds['salt'][:].squeeze()
    th = ds['temp'][:].squeeze()

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
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_splash(in_dict):
    # pretty picture for use in Google Earth

    # START
    fig = plt.figure(figsize=(8,12)) # for real
    #fig = plt.figure(figsize=(8,8)) # for development on laptop
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'temp'#'phytoplankton'
    cmap = 'RdYlBu_r'#'nipy_spectral'
    # things about color limits
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(111)
    # colormap contenders: brg, nipy_spectral, winter, twilight_shifted, ocean_r
    # cm.curl
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=cmap, vlims_fac=2.5)
    pfun.add_coast(ax)
    ax.axis([-129, -122, 42.5, 51.5])
    #ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()

    fig.tight_layout()

    # FINISH
    ds.close()
    plt.show()
    #plt.savefig('/Users/pm7/Desktop/splash.png', transparent=True)

def P_color(in_dict):
    # Code to explore cool color maps

    # START
    fig = plt.figure(figsize=(12,8)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = pfun.get_aa(ds) #[-123.5, -122, 47,48.5]
    ii = 1
    for vn in vn_list:

        # things about color limits
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            vlims_fac=1
        elif vn == 'temp':
            vlims_fac=2.5

        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                aa=aa, vlims_fac=vlims_fac)

        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc=6)
        fig.colorbar(cs, cax=cbaxes, orientation='vertical')

        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
        elif ii == 2:
            ax.set_yticklabels('')
        ii += 1
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_basic_salish(in_dict):

    # START
    fig = plt.figure(figsize=(18,12)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = [-124, -122, 47, 49]
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)

        if vn == 'salt':
            vlims_fac=1
        elif vn == 'temp':
            vlims_fac=2.5

        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                aa=aa, vlims_fac=vlims_fac)
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper right', borderpad=3)
        cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        cb.ax.tick_params(labelsize=fs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels

        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
            pfun.add_windstress_flower(ax, ds, t_scl=.6,
                t_leglen=0.1, center=(.25, .3))
            pfun.add_bathy_contours(ax, ds, depth_levs = [50, 100])
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
                v_scl=30, v_leglen=2, nngrid=80, center=(.1, .1))
            ax.set_yticklabels([])
        ii += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Chl_DO_salish(in_dict):

    # START
    fig = plt.figure(figsize=(18,12)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    aa = [-124, -122, 47, 49]
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if ii == 1:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap='ocean_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
        elif ii == 2:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    slev=0, cmap='rainbow_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Bottom %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
            ax.set_yticklabels([])
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper right', borderpad=3)
        cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        cb.ax.tick_params(labelsize=fs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels

        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)


        ii += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Chl_DO_JdFCanyon(in_dict):

    # START
    fig = plt.figure(figsize=(18,12)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    aa = [-126, -124, 47, 49]
    fs = 14
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if ii == 1:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap='ocean_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
        elif ii == 2:
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    slev=0, cmap='rainbow_r', fac=pinfo.fac_dict[vn])
            ax.set_title('Bottom %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs+2)
            ax.set_yticklabels([])
            pfun.add_bathy_contours(ax, ds, depth_levs = [100, 150, 200, 250, 300, 350, 400], txt=False)
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=10, v_leglen=0.5, nngrid=600, zlev='bot', center=(.8,.05))
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper right', borderpad=3)
        cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        cb.ax.tick_params(labelsize=fs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.tick_params(labelsize=fs) # tick labels

        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)

        ii += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_debug(in_dict):
    # Focused on debugging

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

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
            u = ds['u'][0,-1,:,:].squeeze()
            umax, ujmax, uimax, umin, ujmin, uimin = pfun.maxmin(u)
            v = ds['v'][0,-1,:,:].squeeze()
            vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
            eta = ds['zeta'][0,:,:].squeeze()
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
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_2D(in_dict):
    # For 2D fields.

    # START
    fig = plt.figure(figsize=(20,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['zeta', 'ubar', 'vbar']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict)
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('%s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_dive_vort(in_dict):
    # plots surface fields of divergence and vorticity.

    # currently the color limits are optimized for the sji0 grid, so run as
    # run pan_plot -lt snapshot -pt P_dive_vort -g sj0 -t v0 -x lo8nest -0 2019.08.01

    # START
    fig = plt.figure(figsize=(16,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    # calculate divergence and vorticity
    u = ds['u'][0, -1, :, :]
    v = ds['v'][0, -1, :, :]
    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
    dive = ((np.diff(u, axis=1)/G['DX'][:, 1:-1])[1:-1, :]
            + (np.diff(v, axis = 0)/G['DY'][1:-1, :])[:, 1:-1])
    x = G['lon_psi'] # matrix
    y = G['lat_psi'] # matrix
    dxp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DX'])
    dyp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DY'])
    vort = np.diff(v, axis=1)/dxp - np.diff(u, axis=0)/dyp
    aa = pfun.get_aa(ds)
    scl = 5e-3
    # panel 1
    ax = fig.add_subplot(121)
    cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], dive/scl, cmap='bwr',
                        vmin=-1, vmax=1)
    tstr = ('Surface Divergence (%0.1e $s^{-1}$)' % (scl))
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])
    #
    # panel 2
    ax = fig.add_subplot(122)
    cs = plt.pcolormesh(G['lon_rho'], G['lat_rho'], vort/scl, cmap='bwr',
                        vmin=-1, vmax=1)
    tstr = ('Surface Vorticity (%0.1e $s^{-1}$)' % (scl))
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_willapa_omega2(in_dict):
    # plot bottom and top Omega_arag for Willapa Bay, using PyCO2SYS.

    # START
    fig = plt.figure(figsize=(12,9))
    ds = nc.Dataset(in_dict['fn'])

    verbose=False

    from PyCO2SYS import CO2SYS
    import seawater as sw
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    from time import time

    # specify a geographic region
    aa = [-124.4, -123.6, 46, 47.2] # Willapa Bay

    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
    # find indices that encompass region aa
    i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
    i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
    j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
    j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
    plon = G['lon_psi'][j0:j1-1, i0:i1-1]
    plat = G['lat_psi'][j0:j1-1, i0:i1-1]
    lat = G['lat_rho'][j0:j1,i0:i1] # used in sw.pres

    tt0 = time()
    for ii in [0,1]:
        # loop over bottom and surface layers
        # first extract needed fields and save in v_dict
        v_dict = {}
        vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
        for vn in vn_in_list:
            if ii == 0:
                # depth = bottom
                Ld = G['h'][j0:j1,i0:i1]
                L = zfun.fillit(ds[vn][0,0,j0:j1,i0:i1])
            elif ii==1:
                # depth = surface
                Ld = 0 * G['h'][j0:j1,i0:i1]
                L = zfun.fillit(ds[vn][0,-1,j0:j1,i0:i1])
            v_dict[vn] = L
        # ------------- the CO2SYS steps -------------------------
        # create pressure
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
        PH = CO2dict['pHout']
        PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        ARAG = zfun.fillit(ARAG.reshape((v_dict['salt'].shape)))
        # --------------------------------------------------------

        # PLOT CODE
        fs = 18
        ax = fig.add_subplot(1, 2, ii+1)
        cs = ax.pcolormesh(plon, plat, ARAG[1:-1, 1:-1], vmin=0, vmax=3, cmap='coolwarm_r')
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude', fontsize=fs)
        if ii == 0:
            pfun.add_bathy_contours(ax, ds, depth_levs=[10, 20, 30], txt=True)
            ax.set_title('Bottom Omega',
                fontsize=fs)
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
        elif ii == 1:
            cb = fig.colorbar(cs)
            cb.ax.tick_params(labelsize=fs-2)
            ax.set_title('Surface Omega',
                fontsize=fs)
            ax.set_yticklabels('')
            ax.text(-123.9, 46.85, 'Grays\nHarbor', fontsize=fs-2, style='italic')
            ax.text(-123.86, 46.58, 'Willapa\nBay', fontsize=fs-2, style='italic')
            ax.text(-123.9, 46.04, 'Columbia\nRiver', fontsize=fs-2, style='italic')
        ax.tick_params(labelsize=fs-2)

    if verbose:
        print('   -- carbon took %0.2f sec' % (time()-tt0))

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Carbon(in_dict):

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['alkalinity', 'TIC']
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
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_bio(in_dict):

    # START
    fig = plt.figure(figsize=(20, 8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = vn_list = ['NO3', 'phytoplankton', 'alkalinity', 'TIC']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn in []:
            slev = 0
            ttag = 'Bottom'
        else:
            slev = -1
            ttag = 'Surface'
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        if True:
            ax.axis(pfun.get_aa(ds))
        else:
            aa = [-127.2, -123.8, 45.5, 49.8]
            ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('%s %s %s' % (ttag,pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
        elif ii == 2:
            pfun.add_windstress_flower(ax, ds, center=(.2,.25))
        ii += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer(in_dict):

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    z_level = -30
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
        # cb.formatter.set_useOffset(False)
        # cb.update_ticks()
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
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=1, v_leglen=0.5,
                                      nngrid=80, zlev=z_level)
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer_JdF_Eddy(in_dict):

    # START
    fig = plt.figure(figsize=(24,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa = [-125.8, -124.2, 48, 49]
    pinfo.vlims_dict['salt'] = (32,34)
    pinfo.vlims_dict['temp'] = (6,10)

    vn_list = ['salt', 'temp']
    z_level = -50
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        # if in_dict['auto_vlims']:
        #     pinfo.vlims_dict[vn] = ()
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
        # cb.formatter.set_useOffset(False)
        # cb.update_ticks()
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_bathy_contours(ax, ds, depth_levs = [50, 100, 200], txt=True)
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds, t_scl=0.6, t_leglen=0.1, center=(.15,.75))
        if ii == 2:
            # pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=5, v_leglen=0.2,
            #                           nngrid=30, zlev=z_level)
            pass
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect(in_dict):
    # plots a map and a section (distance, z)

    fs=12
    plt.rc('font', size=fs)

    # START
    fig = plt.figure(figsize=(15,6))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    vn = 'phytoplankton'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    # override
    if vn == 'phytoplankton':
        pinfo.vlims_dict[vn] = (0,25)
        pinfo.vlims_dict['sect_'+vn] = (0,25)
        pinfo.cmap_dict[vn] = 'Spectral_r'

    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -3500
        #x = np.linspace(lon.min(), -124, 500)
        if True:
            x = np.linspace(lon.min(), lon.max(), 500)
            y = 47 * np.ones(x.shape)
        else:
            y = np.linspace(lat.min(), lat.max(), 500)
            x = -126 * np.ones(y.shape)
    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] + 'tracks_new/'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path + track
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
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    #fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    aaf = [-125.5, -122.1, 46.8, 50.3] # focus domain
    #ax.axis(pfun.get_aa(ds))
    ax.axis(aaf)
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    #pfun.add_info(ax, in_dict['fn'])
    #pfun.add_windstress_flower(ax, ds)
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    #
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    # cs = ax.contour(v3['distf'], v3['zrf'], sf,
    #     np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
    #     colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

    plt.rcdefaults()

def P_sect2(in_dict):
    # Plots a map and a section (distance, z) for 2 variables.

    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    #tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
    tracks = ['Line_jdf_v0.p', 'Line_HC_thalweg_long.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
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

    # PLOTTING
    #vn_list = ['NO3', 'phytoplankton']
    vn_list = ['NO3', 'oxygen']
    #vn_list = ['PH', 'ARAG']
    # override colors
    pinfo.vlims_dict['NO3'] = (0,40)
    pinfo.vlims_dict['phytoplankton'] = (0,20)
    pinfo.vlims_dict['oxygen'] = (0,8)
    pinfo.vlims_dict['PH'] = (7, 8)
    pinfo.vlims_dict['ARAG'] = (0, 3)
    counter = 0
    for vn in vn_list:
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

        # map with section line
        if counter == 0:
            ax = fig.add_subplot(2, 3, 1)
        elif counter == 1:
            ax = fig.add_subplot(2, 3, 4)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        pfun.add_coast(ax)
        ax.axis([-126, -122, 47, 50])
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_ylabel('Latitude')
        if counter == 1:
            ax.set_xlabel('Longitude')
        # add section track
        ax.plot(x, y, '-w', linewidth=2)
        ax.plot(x, y, '-k', linewidth=.5)
        ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
            markeredgecolor='r', markeredgewidth=2)
        #
        # section
        if counter == 0:
            ax = fig.add_subplot(2, 3, (2, 3))
        elif counter == 1:
            ax = fig.add_subplot(2, 3, (5, 6))
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(dist.min(), dist.max())
        ax.set_ylim(zdeep, 5)
        sf = pinfo.fac_dict[vn] * v3['sectvarf']
        # set section color limits
        vmin = pinfo.vlims_dict[vn][0]
        vmax = pinfo.vlims_dict[vn][1]
        # plot section
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                           vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
        fig.colorbar(cs)
        cs = ax.contour(v3['distf'], v3['zrf'], sf,
            np.linspace(np.floor(vmin), np.ceil(vmax), 20),
            colors='k', linewidths=0.5)
        if counter == 1:
            ax.set_xlabel('Distance (km)')
            pfun.add_info(ax, in_dict['fn'])
        ax.set_ylabel('Z (m)')
        ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))

        counter += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_willapa(in_dict):
    # Focus on Willapa

    # override
    in_dict['auto_vlims'] = False

    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    vn = 'phytoplankton'
    #vn = 'salt'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    else:
        pinfo.vlims_dict['sect_'+vn] = pinfo.vlims_dict[vn]

    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    track = 'Line_willapa1.p'
    zdeep = -30
    xx = np.array([])
    yy = np.array([])
    track_fn = tracks_path + track
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

    # full map
    ax = fig.add_subplot(1, 4, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)

    # focus map
    ax = fig.add_subplot(1, 4, 2)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    ax.axis([-124.4, -123.6, 46, 47.2])
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 4, (3, 4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_columbia(in_dict):
    # Focus on Columbia River

    # override
    in_dict['auto_vlims'] = False

    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    #vn = 'phytoplankton'
    vn = 'salt'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    else:
        pinfo.vlims_dict[vn] = (0,35)
        pinfo.vlims_dict['sect_'+vn] = (0,35)

    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    track = 'CR_thalweg.p'
    zdeep = -30
    xx = np.array([])
    yy = np.array([])
    track_fn = tracks_path + track
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

    # full map
    ax = fig.add_subplot(1, 4, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)

    aa = [-124.5, -123.3, 45.7, 47.2]
    ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]],
        [aa[2], aa[2], aa[3], aa[3], aa[2]], '-c', linewidth=2)

    # focus map
    ax = fig.add_subplot(1, 4, 2)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 4, (3, 4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sectA(in_dict):
    # plots a map and several sections
    # designed for analytical runs, like aestus1
    # used for MacCready et al. (2018 JPO) variance paper

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOTTING
    vn = 'salt'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    # map and section lines
    ax1 = fig.add_subplot(3,1,1)
    cs = pfun.add_map_field(ax1, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn])
    vlims = pinfo.vlims_dict[vn]
    fig.colorbar(cs)
    ax1.axis([-.5, 1, 44.8, 45.2])
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    pfun.add_info(ax1, in_dict['fn'], fs=9)
    #
    # thalweg section
    x = np.linspace(-0.5, 1, 500)
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    ax = fig.add_subplot(3,1,2)
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-25, 2)
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(1, 35, 35), colors='k', linewidths=.5,)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(vn)
    # add line to map plot
    ax1.plot(x, y, '-k', linewidth=2)
    ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
        markeredgecolor='k', markeredgewidth=2)
    # channel cross-sections
    for ii in range(3):
        y = np.linspace(44.9, 45.1, 200)
        x = 0.3*ii * np.ones(y.shape)
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
        ax = fig.add_subplot(3,3,ii+7)
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(0, 22.5)
        ax.set_ylim(-25, 2)
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
            np.linspace(1, 35, 35), colors='k', linewidths=.5,)
        ax.set_xlabel('Distance (km)')
        if ii==0:
            ax.set_ylabel('Z (m)')
        # add line to map plot
        ax1.plot(x, y, '-k', linewidth=2)
        ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
            markeredgecolor='k', markeredgewidth=2)
    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_basic_sji(in_dict):
    # focus on the San Juan Islands, for the sj0 grid
    # run pan_plot -lt snapshot -pt P_basic_sji -g sj0 -t v0 -x lo8nest -0 2019.08.01

    # START
    fig = plt.figure(figsize=(16,10))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']

    if True:
        aa = pfun.get_aa(ds)
    else:
        #aa = [-123.3, -122.65, 48.3, 48.8]
        aa = [-123.1, -122.8, 48.42, 48.62]

    ii = 1
    for vn in vn_list:
        if False: #in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        else:
            pinfo.vlims_dict['salt'] = (24,31)
            pinfo.vlims_dict['temp'] = (11,15)
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], aa=aa)
        # Inset colorbar
        # from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        # cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper left', borderpad=3)
        # cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        # cb.ax.tick_params(labelsize=fs1-2)
        fig.colorbar(cs, orientation='horizontal')
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds, t_scl=3,
                t_leglen=0.1, center=(.33, .52))
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
                v_scl=80, v_leglen=2, nngrid=60, center=(.27, .45))
            ax.set_yticklabels([])
        ii += 1
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_speed_sji(in_dict):
    # focus on the San Juan Islands, for the sj0 grid
    # with a cool speed color and some vectors

    # run pan_plot -lt snapshot -pt P_speed_sji -g sj0 -t v0 -x lo8nest -0 2019.08.01

    # START
    fig = plt.figure(figsize=(16,10))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    # get surface velocity
    uu = ds['u'][0, -1, :, :].squeeze()
    vv = ds['v'][0, -1, :, :].squeeze()
    # interpolate to trimmed rho grid
    u = (uu[1:, :] + uu[:-1, :])/2
    v = (vv[:, 1:] + vv[:, :-1])/2
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    speed = np.sqrt(u**2 + v**2)
    # encompasing axes
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]

    aa = pfun.get_aa(ds)

    vn = 'salt'
    #pinfo.vlims_dict['salt'] = (26,31)
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(121)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], aa=aa)
    fig.colorbar(cs, orientation='horizontal')
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds, t_scl=3,
        t_leglen=0.1, center=(.33, .52))

    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x, y, speed, vmin=0, vmax=2, cmap='rainbow')
    fig.colorbar(cs, orientation='horizontal')
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface speed (m/s)')
    ax.set_xlabel('Longitude')
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
        v_scl=80, v_leglen=2, nngrid=60, center=(.27, .45))
    ax.set_yticklabels([])

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_sji0(in_dict):
    # Plots a map and a section (distance, z) for 2 variables.

    # run pan_plot -lt snapshot -pt P_sect_sji0 -g sj0 -t v0 -x lo8nest -0 2019.08.01

    # START
    fig = plt.figure(figsize=(18,10))
    ds = nc.Dataset(in_dict['fn'])
    
    from cmocean import cm
    

    # PLOT CODE
    #
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    zdeep = -200
    tracks_path = Ldir['data'] + 'tracks_new/'
    track_fn = tracks_path + 'sji_thalweg0.p'
    # get the track to interpolate onto
    pdict = pickle.load(open(track_fn, 'rb'))
    xx = pdict['lon_poly']
    yy = pdict['lat_poly']
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
            if (x0 != x1) or (y0 != y1):
                # only add points if they are not a repeat
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
            else:
                pass

    # PLOTTING
    vn_list = ['salt', 'temp']
    # override colors
    pinfo.vlims_dict['salt'] = (28,32)
    pinfo.vlims_dict['temp'] = (8,16)
    ci = 0.2 # contour interval
    pinfo.cmap_dict['salt'] = cm.haline
    pinfo.cmap_dict['temp'] = cm.thermal
    aa = pfun.get_aa(ds)
    counter = 0
    for vn in vn_list:
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
        # find indices of some sills
        sind_list = []
        sdist_list = [6.27,25.82,39.28]
        smarker_list = ['o','*','^']
        for sdist in sdist_list:
            sind_list.append(zfun.find_nearest_ind(dist,sdist))
        # map with section line
        if counter == 0:
            ax = fig.add_subplot(2, 3, 1)
        elif counter == 1:
            ax = fig.add_subplot(2, 3, 4)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_ylabel('Latitude')
        if counter == 1:
            ax.set_xlabel('Longitude')
        # add section track
        ax.plot(x, y, '-w', linewidth=2)
        ax.plot(x, y, '-k', linewidth=.5)
        # ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        #     markeredgecolor='r', markeredgewidth=2)
        scounter = 0
        for sind in sind_list:
            ax.plot(x[sind], y[sind], marker=smarker_list[scounter],
                markersize=7, markerfacecolor='w',
                markeredgecolor='k', markeredgewidth=1)
            scounter += 1
        #
        # section
        if counter == 0:
            ax = fig.add_subplot(2, 3, (2, 3))
        elif counter == 1:
            ax = fig.add_subplot(2, 3, (5, 6))
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(dist.min(), dist.max())
        ax.set_ylim(zdeep, 5)
        sf = pinfo.fac_dict[vn] * v3['sectvarf']
        # set section color limits
        vmin = pinfo.vlims_dict[vn][0]
        vmax = pinfo.vlims_dict[vn][1]
        # plot section
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                           vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
        fig.colorbar(cs)
        cs = ax.contour(v3['distf'], v3['zrf'], sf,
            np.arange(np.floor(vmin), np.ceil(vmax), .2),
            colors='k', linewidths=0.5)
        scounter = 0
        for sind in sind_list:
            ax.plot(dist[sind], 0, marker=smarker_list[scounter],
                markersize=10, markerfacecolor='w',
                markeredgecolor='k', markeredgewidth=1)
            scounter += 1
        # ax.plot(dist[is0], 0, '*r', markersize=5, markerfacecolor='w',
        #     markeredgecolor='r', markeredgewidth=2)
        if counter == 1:
            ax.set_xlabel('Distance (km)')
            pfun.add_info(ax, in_dict['fn'])
        ax.set_ylabel('Z (m)')
        ax.set_title('Section %s %s (C.I.=%0.1f)' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn],ci))

        counter += 1

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
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
    fn = in_dict['fn']
    ds = nc.Dataset(fn)
    
    gtagex = in_dict['fn'].split('/')[-3]
    year_str = in_dict['fn'].split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] + 'superplot/forcing_'+gtagex+'_'+year_str+'.p'
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] + 'tracks_new/'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
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
    lon = ds['lon_psi'][:]
    lat = ds['lat_psi'][:]
    v =ds[vn][0, -1, 1:-1, 1:-1]
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
        + datetime.strftime(T['tm'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
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
    tm = T['tm'] # datetime
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
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_oxygen(in_dict):
    # Plot bottom oxygen maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'oxygen'
    vlims = (0, 8) # full map
    vlims2 = (0, 8) # PS map
    vlims3 = (0, 8) # PS section
    from cmocean import cm
    cmap = cm.oxy

    # get model fields
    fn = in_dict['fn']
    ds = nc.Dataset(fn)
    
    gtagex = in_dict['fn'].split('/')[-3]
    year_str = in_dict['fn'].split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] + 'superplot/forcing_'+gtagex+'_'+year_str+'.p'
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] + 'tracks_new/'
    tracks = ['Line_HC_thalweg_long.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
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
    lon = ds['lon_psi'][:]
    lat = ds['lat_psi'][:]
    v =ds[vn][0, 0, 1:-1, 1:-1]
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
    ax.text(.95, .03, 'LiveOcean\nBottom Oxygen\n'
        + datetime.strftime(T['tm'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    # ax.text(.99,.97,'DO range\n'+ str(vlims), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=.85*fs)
    ax.text(1, .85, r'$[mg\ L^{-1}]$', transform=ax.transAxes, fontsize=fs, ha='right')

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
    tm = T['tm'] # datetime
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
    if len(in_dict['fn_out']) > 0:
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
    fn = in_dict['fn']
    ds = nc.Dataset(fn)
    
    gtagex = in_dict['fn'].split('/')[-3]
    year_str = in_dict['fn'].split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] + 'superplot/forcing_'+gtagex+'_'+year_str+'.p'
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] + 'tracks_new/'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
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
    lon = ds['lon_psi'][:]
    lat = ds['lat_psi'][:]
    v =ds[vn][0, -1, 1:-1, 1:-1]
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, -6:] = np.nan
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
        + datetime.strftime(T['tm'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.99,.97,'range\n'+ str(vlims), transform=ax.transAxes,
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
    # ax.text(.93,.97,'range\n'+ str(vlims2), transform=ax.transAxes,
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
    # ax.text(.99,.4,'range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['tm'] # datetime
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
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
