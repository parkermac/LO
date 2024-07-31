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
import sys

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

from cmocean import cm # have to import after matplotlib to work on remote machine
    
def P_basic(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
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
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
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
        
def P_basic_AttSW_boxes(in_dict):
    # One-off code to look at the places in Aurora's code where AttSW is increased
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
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
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
            
            # AttSW boxes
            pfun.draw_box(ax, [-125.31, -123.89, 49.13, 51.02])
            pfun.draw_box(ax, [-123.89, -122, 47.02, 50.29])
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_basic_PS(in_dict):
    # Like P_basic but focused on Puget Sound
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = [-123.2, -122.1, 47, 48.4]
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
    pinfo.cmap_dict['salt'] = 'jet'
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            # this gets the color limits to be focused on the "aa" region
            i0 = zfun.find_nearest_ind(ds.lon_rho[0,:],aa[0])
            i1 = zfun.find_nearest_ind(ds.lon_rho[0,:],aa[1])
            j0 = zfun.find_nearest_ind(ds.lat_rho[:,0],aa[2])
            j1 = zfun.find_nearest_ind(ds.lat_rho[:,0],aa[3])
            fld = ds[vn][0,-1,j0:j1,i0:i1].values
            vlims = pfun.auto_lims(fld)
            pinfo.vlims_dict[vn] = vlims
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Fb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'Fb']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'salt':
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                    vlims_fac=pinfo.range_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
                    fontsize=1.2*fs)
        elif vn == 'Fb':
            # plot vertically integrateed buoyancy flux
            # calculate potential density
            import seawater as sw
            rho = sw.dens0(ds.salt.values.squeeze(), ds.temp.values.squeeze())
            # calculate vertically-integrated buoyancy flux
            Fb = -9.8 * np.sum(ds.AKs[0,1:-1,:,:].squeeze() * np.diff(rho, axis=0), axis=0).values
            Fb[ds.mask_rho.values==0] = np.nan
            plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
            cs = ax.pcolormesh(plon,plat, Fb, vmin=0, vmax=.5, cmap='YlOrRd')
            ax.set_title('Vertically Integrated Fb [W/m2]')
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        elif ii == 2:
            ax.set_yticklabels([])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_fancy(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    # PLOT CODE
    import cmcrameri.cm as cmr
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            cmap = cmr.hawaii_r
            vlims_fac = .5
        elif vn == 'temp':
            cmap = cmr.roma_r
            vlims_fac = 1
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=cmap, fac=pinfo.fac_dict[vn], vlims_fac=vlims_fac)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 8.5
    pfun.start_plot(fs=fs, figsize=(hgt*2.7/AR,hgt))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 2*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(1, 2, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            pass
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort2(in_dict):
    # same as dive_vort but focused on a specific region
    
    # JdF:
    aa = [-125, -122.3, 47.8, 48.8]
    
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    # aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 6
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 4*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(2, 1, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_ylabel('Latitude')
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            #pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_xlabel('Longitude')
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
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
    pfun.start_plot(fs=fs, figsize=(18,10))
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
        
def P_bot_top_arag(in_dict):
    # START
    aa = [-124.4,-123.7,46.35,47.1]  # hard-coded for Willapa & Grays
    fs = 14
    pfun.start_plot(fs=fs, figsize=(18,10))
    fig = plt.figure()
    fn = in_dict['fn']
    ds = xr.open_dataset(fn)
    arag_bot, arag_top, px, py = pfun.get_bot_top_arag(fn, aa=aa)
    # PLOT CODE
    fs = 14
    for ii in [1,2]:
        if ii == 1:
            fld = arag_bot
            stext = 'Bottom'
        elif ii == 2:
            fld = arag_top
            stext = 'Surface'
        ax = fig.add_subplot(1, 2, ii)
        cs = ax.pcolormesh(px,py,fld, vmin=0, vmax=3, cmap='coolwarm_r')
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title(r'%s $\Omega_{arag}$' % (stext), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        pfun.add_bathy_contours(ax, ds, txt=True)
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_DO_WA_shelf(in_dict):
    # Focus on bottom DO on the WA shelf
    aa = [-126.1, -123.7, 45.8, 48.8]
    xtl = [-126, -125, -124]
    ytl = [46, 47, 48]
    
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(7,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'oxygen'
    slev = 0
    stext = 'Bottom'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(111)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
    fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
    ax.set_xlabel('Longitude')
    pfun.add_bathy_contours(ax, ds, txt=False)
    ax.set_ylabel('Latitude')
    ax.set_xticks(xtl)
    ax.set_yticks(ytl)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    
    pfun.add_windstress_flower(ax, ds, t_scl=0.5, t_leglen=0.1, center=(.85,.65), fs=12)
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
        
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
    vn_list = ['u', 'v', 'zeta']
    do_wetdry = False
    
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(8*len(vn_list),10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    ii = 1
    for vn in vn_list:
        if 'lon_rho' in ds[vn].coords:
            tag = 'rho'
        if 'lon_u' in ds[vn].coords:
            tag = 'u'
        if 'lon_v' in ds[vn].coords:
            tag = 'v'
        x = ds['lon_'+tag].values
        y = ds['lat_'+tag].values
        px, py = pfun.get_plon_plat(x,y)
        if vn in ['u', 'v']:
            v = ds[vn][0,-1,:,:].values
            vmin = -2
            vmax = 2
            cmap='hsv_r'
        elif vn == 'zeta':
            v = ds[vn][0,:,:].values
            h = ds.h.values
            mr = ds.mask_rho.values
            v[mr==0] = np.nan
            h[mr==0] = np.nan
            v = v + h
            vn = 'depth'
            vmin = 2
            vmax = 4
            cmap='RdYlGn'
        else:
            v = ds[vn][0, -1,:,:].values
        ax = fig.add_subplot(1, len(vn_list), ii)
        ax.set_xticks([])
        ax.set_yticks([])
        cs = ax.pcolormesh(px, py, v, cmap=cmap, vmin=vmin, vmax=vmax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'], his_num=True)
        vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
        ax.plot(x[vjmax,vimax], y[vjmax,vimax],'*y', mec='k', markersize=15)
        ax.plot(x[vjmin,vimin], y[vjmin,vimin],'oy', mec='k', markersize=10)
        ax.set_title(('%s ((*)max=%0.1f, (o)min=%0.1f)' % (vn, vmax, vmin)))
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
    #vn_list = ['oxygen', 'temp']
    vn_list = ['TIC', 'alkalinity']
    z_level = -1500
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

def P_layer_CUC(in_dict):
    """
    Focused on plotting a layer relevant to seeing the California Undercurrent.

    This is also the first use of pfun.gat_laym_alt() which uses np.argmax()
    and np.take_along_axis().
    
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,8))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    from cmocean import cm
    # from time import time
    vn_list = ['salt', 'speed']
    # note that 'speed' is a derived quantity
    z_level = -400
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'speed':
            zfull_u = pfun.get_zfull(ds, in_dict['fn'], 'u')
            zfull_v = pfun.get_zfull(ds, in_dict['fn'], 'v')
            # tt0 = time()
            u = pfun.get_laym_alt(ds, zfull_u, ds.mask_u.values, 'u', z_level)
            v = pfun.get_laym_alt(ds, zfull_v, ds.mask_v.values, 'v', z_level)
            # print('time for u,v layers = %0.2f sec' % (time()-tt0))
            # interpolate to rho grid
            ur = np.nan * np.ones(ds.mask_rho.shape)
            vr = np.nan * np.ones(ds.mask_rho.shape)
            ur[1:-1,1:-1] = u[1:-1,:-1] + (u[1:-1,1:]-u[1:-1,:-1])/2
            vr[1:-1,1:-1] = v[:-1,1:-1] + (v[1:,1:-1]-v[:-1,1:-1])/2
            # create speed and mask
            speed = np.sqrt(ur*ur + vr*vr)
            speed[ds.mask_rho.values == 0] = np.nan
            v_scaled = speed
            speed_max = 0.2
            vlims = (0,speed_max)
            pinfo.cmap_dict['speed'] = cm.amp #cm.speed
            pinfo.tstr_dict['speed'] = 'Speed'
            pinfo.units_dict['speed'] = ' $(m\ s^{-1})$'
        else:
            laym = pfun.get_laym_alt(ds, zfull, ds.mask_rho.values, vn, z_level)
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
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level, v_leglen=speed_max, v_scl=5)
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

def P_bpress(in_dict):
    """
    Specialized plot related to bottom pressure anomalies.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    sta_dict = {
        'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
        'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
        'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
        'PN01A':(-125.3983, 44.5096), # Slope Base (2905 m)
        }
    vn_list = ['salt', 'temp']
    z_level = -300
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
            if vn == 'salt':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=0.3)
            elif vn == 'temp':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=2)
            else:
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
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level, v_scl=5, v_leglen=0.3)
            for sta in sta_dict.keys():
                ax.plot(sta_dict[sta][0], sta_dict[sta][1], 'o', mfc='y', mec='k', ms=10)
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
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'#'phytoplankton'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'jet'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -300
        x_e = np.linspace(-124, -123, 500)
        y_e = 48.368 * np.ones(x_e.shape)
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
                x_e = np.linspace(x0, x1, nn)
                y_e = np.linspace(y0,y1, nn)
            else:
                x_e = np.concatenate((x_e, np.linspace(x0, x1, nn)[1:]))
                y_e = np.concatenate((y_e, np.linspace(y0, y1, nn)[1:]))
                
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)
        
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * fld_s
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
    aaf = [-123.5, -122.3, 48, 49.5] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    
    ax.plot(dist_se[0,:], zw_se[0,:], '-k', linewidth=2)
    ax.plot(dist_se[-1,:], zw_se[-1,:], '-k', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(dist_se,zw_se,sf,
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
        
def P_sect_CR(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    Uses the new pfun.get_sect() function.

    This version is focused on the Columbia River Estuary.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,8))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'#'phytoplankton'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'Spectral_r'#cm.haline
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['CR_thalweg.p']
    zdeep = -30
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
            x_e = np.linspace(x0, x1, nn)
            y_e = np.linspace(y0,y1, nn)
        else:
            x_e = np.concatenate((x_e, np.linspace(x0, x1, nn)[1:]))
            y_e = np.concatenate((y_e, np.linspace(y0, y1, nn)[1:]))
                
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # Convert dist to be in units of "River Mile" (RM) to enable
    # easier comparinson to the figures in Jay and Smith (1990).
    # We assume that RM = 0 between the jetties, which is at dist = 24 km
    # for the section we use here.
    km0 = 24 # km into dist to get to RM 0.
    km2mi = 0.621371 # convert km to miles
    dist = (dist - km0)*km2mi
    dist_e = (dist_e - km0)*km2mi
    dist_se = (dist_se - km0)*km2mi
        
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * fld_s
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
    aaf = [-125, -123.5, 45, 47] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    # add a marker on the section near the mouth
    nm = np.argmax(dist>0)
    ax.plot(x[nm], y[nm], 'or', markersize=5, markerfacecolor='k',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-125, -124])
    ax.set_yticks([45, 46, 47])

    # pfun.add_bathy_contours(ax, ds, depth_levs = [5,7.5,10], txt=True)

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    
    ax.plot(dist_e, zw_se[0,:], '-k', linewidth=2)
    ax.plot(dist_e, zw_se[-1,:], '-k', linewidth=1)
    ax.set_xlim(dist_e.min(), dist_e.max())
    ax.set_ylim(zdeep, 5)

    # add a marker on the section near the mouth
    nm = np.argmax(dist_e>0)
    ax.plot(dist_e[nm], zw_se[-1,nm], 'or', markersize=5, markerfacecolor='k',
        markeredgecolor='r', markeredgewidth=2)
    ax.plot(dist_e[0], zw_se[-1,0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(dist_se,zw_se,sf,
        vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    # add contours
    # note: we have to make x and z arrays of the same shape as sf
    dist_s = (dist_se[:,0:-1] + np.diff(dist_se,axis=1)/2)[:-1,:]
    zw_s = zw_se[:,0:-1] + np.diff(zw_se,axis=1)/2
    zr_s = zw_s[0:-1,:] + np.diff(zw_s,axis=0)/2
    levs = np.arange(5,35,5)
    ax.contour(dist_s,zr_s,sf,levs,colors='w',linewidths=2)
    ax.text(.95,.05,'Salinity contours: %s' % (str(levs)),ha='right',
        transform=ax.transAxes)

    ax.set_xlabel('Distance (River Mile)')
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
        
def P_sect_hc(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    Focus on Hood Canal.
    
    Uses the new pfun.get_sect() function.
    
    2023.10.31 This is the first code to use the new location of user-created sections
    in LO_output/section_lines, created using the new tool LO/plotting/create_sect_lines.py.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'oxygen'#'phytoplankton'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'jet'
    elif vn == 'oxygen':
        # pinfo.cmap_dict[vn] = 'Spectral_r'
        pinfo.cmap_dict[vn] = 'nipy_spectral'
        # pinfo.cmap_dict[vn] = 'tab20b'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # read in a section (or list of sections)
    tracks_path = Ldir['LOo'] / 'section_lines'
    tracks = ['hc1.p','hc2.p','hc3.p']
    zdeep = -175
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        tdf = pd.read_pickle(track_fn)
        xx = np.concatenate((xx,tdf['x'].to_numpy()))
        yy = np.concatenate((yy,tdf['y'].to_numpy()))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x_e = np.linspace(x0, x1, nn)
            y_e = np.linspace(y0,y1, nn)
        else:
            x_e = np.concatenate((x_e, np.linspace(x0, x1, nn)[1:]))
            y_e = np.concatenate((y_e, np.linspace(y0, y1, nn)[1:]))
                
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)
        
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * fld_s
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
    aaf = [-123.5, -122.5, 47.25, 48.5] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-k', linewidth=2)
    ax.plot(x[0], y[0], 'o', markersize=5, markerfacecolor='w',
        markeredgecolor='k', markeredgewidth=2)
    ax.set_xticks([-123.5, -123, -122.5])
    ax.set_yticks([47.5, 48, 48.5])

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    
    ax.plot(dist_se[0,:], zw_se[0,:], '-k', linewidth=2)
    ax.plot(dist_se[-1,:], zw_se[-1,:], '-k', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(dist_se,zw_se,sf,
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
        
def P_sect_soundspeed(in_dict):
    """
    Soundspeed section plot
    """
    import gsw

    ds = xr.open_dataset(in_dict['fn'])
    # create track by hand
    x_e = np.linspace(-124.85,-124.2, 100) # shelf only
    #x = np.linspace(-126,-124.2, 100) # shows SOFAR channel
    y_e = 47 * np.ones(x_e.shape)
    
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, s_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], 'salt', x_e, y_e)
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, th_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], 'temp', x_e, y_e)
    
    # make Z on section centers, with surface at 0
    Z = zw_se - zw_se[-1,:]
    Z = Z[:-1,:] + np.diff(Z, axis=0)/2
    Z = Z[:,:-1] + np.diff(Z, axis=1)/2

    p = gsw.p_from_z(Z, 47)
    SA = gsw.SA_from_SP(s_s, p, -125, 47)
    CT = gsw.CT_from_pt(SA, th_s)
    spd = gsw.sound_speed(SA, CT, p)

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(16,9))
    fig, axes = plt.subplots(nrows=3, ncols=2)

    ax = axes[0,0]
    cs = ax.pcolormesh(dist_se, zw_se, SA, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Absolute Salinity', transform=ax.transAxes, ha='right')

    ax = axes[1,0]
    cs = ax.pcolormesh(dist_se, zw_se, CT, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,0]
    cs = ax.pcolormesh(dist_se, zw_se, spd, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    ax = axes[0,1]
    ax.plot(SA,Z, alpha=.2)
    ax.text(.05, .05, 'Absolute Salinity', transform=ax.transAxes, ha='left')

    ax = axes[1,1]
    ax.plot(CT,Z, alpha=.2)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,1]
    ax.plot(spd,Z, alpha=.2)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    fig.suptitle(str(in_dict['fn']))

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
    import gsw
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
        lon = G['lon_rho'][j0:j1,i0:i1] # used in gsw.SA_from_SP
        lat = G['lat_rho'][j0:j1,i0:i1] # used in gsw.SA_from_SP
        # first extract needed fields and save in v_dict
        v_dict = {}
        # This version calculates in situ density using the gsw toolbox.
        # The results were identidacal to those using the roms rho and sw to within
        # 4 decimal places in ARAG.
        vn_in_list = ['temp', 'salt', 'alkalinity', 'TIC']
        for cvn in vn_in_list:
            L = ds[cvn][0,nlev,j0:j1,i0:i1].values
            v_dict[cvn] = L
        # ------------- the CO2SYS steps -------------------------
        Ld = G['h'][j0:j1,i0:i1]
        if nlev == 0:
            pass
        elif nlev == -1:
            Ld = 0 * Ld
        else:
            print('get_arag() error: nlev must be 0 or -1')
            sys.exit()
        Lpres = gsw.p_from_z(-Ld, lat) # pressure [dbar]
        SA = gsw.SA_from_SP(v_dict['salt'], Lpres, lon, lat)
        CT = gsw.CT_from_pt(SA, v_dict['temp'])
        rho = gsw.rho(SA, CT, Lpres) # in situ density
        Ltemp = gsw.t_from_CT(SA, CT, Lpres) # in situ temperature
        # convert from umol/L to umol/kg using in situ dentity
        Lalkalinity = 1000 * v_dict['alkalinity'] / rho
        Lalkalinity[Lalkalinity < 100] = np.nan
        LTIC = 1000 * v_dict['TIC'] / rho
        LTIC[LTIC < 100] = np.nan
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # PH = CO2dict['pHout']
        # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        ARAG = ARAG.reshape((v_dict['salt'].shape))
        ARAG = ARAG[1:-1, 1:-1]
        # print(np.nanmax(ARAG))
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
    tstr = T['dt'].strftime(Lfun.ds_fmt)
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
        
def P_splash2(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Salish Sea.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values
    ox = pinfo.fac_dict['oxygen'] * ds['oxygen'][0,0,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
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
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Salish Sea
    aa = [-125.3, -122.1, 46.8, 50.3]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    from cmocean import cm
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,ox, cmap=cm.oxy, vmin=0, vmax=10)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    ax.set_xlabel('Longitude')
    ax.text(.84,.95,'Salish Sea\n\nBottom Oxygen\n$[mg\ L^{-1}]$',
         ha='right', va='top', weight='bold', transform=ax.transAxes)
         
    # add labels
    ax.text(-122.8,49.335,'Fraser\nRiver',size=fs2,
        style='italic',ha='center',va='center',rotation=0)
    ax.text(-123.7,49.2528,'Strait of Georgia',size=fs2,
        style='italic',ha='center',va='center',rotation=-30)
    ax.text(-123.5,48.28,'Strait of Juan de Fuca',size=fs2,
        style='italic',ha='center',va='center',rotation=0,
        color='w')
    ax.text(-123.3,47.6143,'Puget\nSound',size=fs2,
        style='italic',ha='center',va='center',rotation=+55)
    ax.text(-122.3,48.48,'Skagit\nRiver',size=fs3,
        style='italic',ha='center',va='center',
        bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(-123.173,48.44,'Haro\nStrait',size=fs3,
        style='italic',ha='center',va='center',
        color='w')

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash3(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    th = ds['temp'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap='gist_earth', shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash4(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound Salinity.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True
    cmap = 'Spectral_r'
    vlims = (25, 33) # full map
    vlims2 = (25, 33) # PS map

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    salt = ds['salt'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 1
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface Salinity $[g\ kg^{-1}]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims2[0], vmax=vlims2[1])
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash_sneaker(in_dict):
    """
    This makes a fancy plot for a Nat Geo request, including some
    drifter tracks with 0.03 windage from 2019.01.18.
    
    Run as:
    run pan_plot.py -0 2019.01.18 -gtx cas6_v0_live -pt P_splash_sneaker
    """
    
    no_text = True # a flag to turn off all text additions, for the Japanese version
    
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True
    cmap = 'RdYlBu_r'
    vn = 'temp'
    vlims = (6,14)

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    salt = ds[vn][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15plus' / 'topo.nc')
        tds = xr.open_dataset(tfn)
        step = 1
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan
        
    # get drifter tracks
    dfn = Ldir['LOo'] / 'tracks' / 'sneaker_surf_wind3' / 'release_2019.01.18.nc'
    dds = xr.open_dataset(dfn)
    skp = 1
    dx = dds.lon[:,::skp].values
    dy = dds.lat[:,::skp].values

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    if no_text:
        pass
    else:
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        tstr = T['dt'].strftime(Lfun.ds_fmt)
        ax.text(.98,.99,'Initial release\nlocations', size=fs*1.5,
             ha='right', va='top', weight='bold', transform=ax.transAxes,
             bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.text(.17,.05,'LiveOcean\nSurface water\nTemperature $[^{\circ}C]$\n\n'+tstr,
             ha='left', va='bottom', weight='bold', transform=ax.transAxes,
             bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for drifter release
    aa = [-125.8, -122.2, 46.3, 51]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='b', alpha=1, linewidth=3, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    if no_text:
        pass
    else:
        ax.text(-123.072,46.7866,'Washington', size=fs2,
            style='italic',ha='center',va='center',rotation=-45)
        ax.text(-122.996,44.5788,'Oregon', size=fs2,
            style='italic',ha='center',va='center',rotation=-45)
    
        ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
            style='italic',ha='center',va='center',rotation=-45)
        ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
            style='italic',ha='center',va='center',rotation=-10)
        
    # add drifter tracks
    ax.plot(dx[0,:],dy[0,:],'.k', ms=1, alpha=.5)
    
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    if no_text:
        pass
    else:
        ax.set_xlabel('Longitude')
        ax.text(.98,.99,'Locations after\nthree days', size=fs*1.5, c='k',
             ha='right', va='top', weight='bold', transform=ax.transAxes,
             bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
         
    # add drifter tracks
    # ax.plot(dx, dy, '-k', lw=.5, alpha=.5)
    ax.plot(dx[-1,:],dy[-1,:],'.k')
    
    plt.setp(ax.spines.values(), linewidth=3, color='b')
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_admiralty(in_dict):
    """
    Designed to explore the TEF salt flux at Admiralty Inlet.
    
    This is quite task-specific code, like the "superplots" and
    only works if you have a tef section already extracted.
    """
    # START
    ds = xr.open_dataset(in_dict['fn'])
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,12))
    fig = plt.figure()
    
    # PLOT CODE
    
    # load a TEF time series
    import sys
    pth = str(Ldir['LO'] / 'extract' /'tef')
    if pth not in sys.path:
        sys.path.append(pth)
    import flux_fun
    import tef_fun
    year_str = in_dict['fn'].parent.name.split('.')[0][1:]
    gtagex = 'cas6_v3_lo8b' #in_dict['fn'].parent.parent.name
    dates_str = '_' + year_str + '.01.01_' + year_str + '.12.31'
    tef_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef' / ('bulk' + dates_str)
    gridname = gtagex.split('_')[0]
    sect_name = 'ai4'
    df, in_sign, _, _ = flux_fun.get_two_layer(tef_dir, sect_name, gridname)
    sect_df = tef_fun.get_sect_df(gridname)
    xs0, xs1, ys0, ys1 = sect_df.loc[sect_name,:]
    
    Qin = df['Qin']/1000
    Qout = df['Qout']/1000
    Qprism = (df['qabs']/2)/1000
    Sin = df['salt_in']
    Sout = df['salt_out']
    DS = Sin - Sout
    tef_day = df.index.dayofyear
    
    # get this time
    ot = ds.ocean_time.values[0]
    tpd = pd.Timestamp(ot)
    dd = tpd.dayofyear + tpd.hour/24
    idd = zfun.find_nearest_ind(tef_day,dd)
    
    sin = Sin[idd]
    sout = Sout[idd]
    
    # plot TEF things
    ax = fig.add_subplot(413)
    ax.plot(tef_day, Sin, '-r',tef_day, Sout, '-b', linewidth=3)
    ax.plot(dd, sin, 'or',dd, sout, 'ob', markersize=10)
    ax.set_xlim(0,365)
    ax.set_xticklabels([])
    ax.axvline(dd)
    ax.text(.05,.8,'$S_{in}$', c='r',fontweight='bold',
        transform=ax.transAxes, fontsize= 24,
        bbox=dict(facecolor='w', edgecolor='None',alpha=.9))
    ax.text(.05,.1,'$S_{out}$', c='b',fontweight='bold',
        transform=ax.transAxes, fontsize= 24,
        bbox=dict(facecolor='w', edgecolor='None',alpha=.9))
    
    ax = fig.add_subplot(414)
    ax.plot(tef_day, Qprism, '-g', linewidth=3)
    ax.set_xlim(0,365)
    ax.set_xlabel('Yearday ' + year_str)
    ax.axvline(dd)
    ax.text(.05,.8,'$Q_{prism}\ [10^{3}m^{3}s^{-1}]$', c='g',fontweight='bold',
        transform=ax.transAxes, fontsize= 24,
        bbox=dict(facecolor='w', edgecolor='None',alpha=.9))
    
    # Thickness maps
    
    # to speed things up we will clip the region of interest
    aa = [-123, -122.2, 47.7, 48.3]
    x0, x1, y0, y1 = aa
    Lon = ds.lon_rho[0,:].values
    Lat = ds.lat_rho[:,0].values
    i0 = zfun.find_nearest_ind(Lon, x0)
    i1 = zfun.find_nearest_ind(Lon, x1)
    j0 = zfun.find_nearest_ind(Lat, y0)
    j1 = zfun.find_nearest_ind(Lat, y1)
    lon = ds.lon_rho[j0:j1,i0:i1].values
    lat = ds.lat_rho[j0:j1,i0:i1].values
    h = ds.h[j0:j1,i0:i1].values
    mask = ds.mask_rho[j0:j1,i0:i1].values
    xp, yp = pfun.get_plon_plat(lon,lat)
    
    salt = ds.salt[0,:,j0:j1,i0:i1].values
    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    zw = zrfun.get_z(h, 0*h, S, only_w=True)
    dz = np.diff(zw, axis=0)
    
    # plot thickness of water with salinity greater than some value
    salt_mask = salt > sin
    dzc = dz.copy()
    dzc[~salt_mask] = 0
    thick = dzc.sum(axis=0)
    thick[mask==0] = np.nan
    ax = fig.add_subplot(221)
    cs = ax.pcolormesh(xp, yp, thick, vmin=0, vmax=50, cmap='PuRd')
    fig.colorbar(cs, ax=ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.text(.95,.85,'Thickness [m]\n$S >\ S_{in}$',
        transform=ax.transAxes, ha='right',
        bbox=dict(facecolor='w', edgecolor='None',alpha=.9))
    # add section
    ax.plot([xs0, xs1], [ys0,ys1],'-g', linewidth=3)
    # add info
    pfun.add_info(ax, in_dict['fn'])
    
    # plot thickness of water with salinity less than some value
    salt_mask = salt < sout
    dzc = dz.copy()
    dzc[~salt_mask] = 0
    thick = dzc.sum(axis=0)
    thick[mask==0] = np.nan
    ax = fig.add_subplot(222)
    cs = ax.pcolormesh(xp, yp, thick, vmin=0, vmax=50, cmap='YlGnBu')
    fig.colorbar(cs, ax=ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.text(.95,.85,'Thickness [m]\n$S <\ S_{out}$',
        transform=ax.transAxes, ha='right',
        bbox=dict(facecolor='w', edgecolor='None',alpha=.9))
    # add section
    ax.plot([xs0, xs1], [ys0,ys1],'-g', linewidth=3)
    
    # fig.tight_layout()
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
    x_e = x
    y_e = y
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)
    
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
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
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
    ax.plot(dist, ztop+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], ztop[n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], ztop[n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * fld_s
    # plot section
    cs = ax.pcolormesh(dist_se, zw_se, sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
        va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
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
    x_e = x
    y_e = y
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

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
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
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
    ax.plot(dist, ztop+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], ztop[n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], ztop[n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * fld_s
    # plot section
    cs = ax.pcolormesh(dist_se, zw_se, sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # labels
    ax.text(0, 0, 'SECTION\nHood Canal', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
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
        
def P_superplot_oxygen_coast(in_dict):
    # Plot bottom oxygen maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.
    
    # 2023.05.19 Based on P_superplot_oxygen but focused on the coast
    
    # Here is the command I used to run it on parigee:
    # run pan_plot.py -gtx cas6_v0_live -ro 0 -0 2019.01.01 -1 2019.12.31 -lt daily -pt P_superplot_oxygen_coast -mov True
    # This takes awhile and would have been better done from the bash command line in background.
    # I had to copy over a superplot mooring extraction for this gtagex from my laptop.

    vn = 'oxygen'
    vlims = (0, 10) # full map
    vlims2 = (0, 10) # coast map
    vlims3 = (0, 10) # coast section
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
    xx = np.array([-125.125, -124.125])
    yy = np.array([47, 47])
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 200
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    x_e = x
    y_e = y
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

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
    aa = [-125.8, -123.5, 46, 49]
    # aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nBottom Oxygen\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=.85*fs)
    ax.text(1, .85, r'$[mg\ L^{-1}]$', transform=ax.transAxes, fontsize=fs, ha='right')

    # coast map
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
    # ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
    #     markeredgecolor='k')
    # ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
    #     markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(233)
    ax.plot(dist, ztop+5, linestyle='--', color='k', linewidth=2)
    # ax.plot(dist[n_ai], ztop[n_ai] + 5, marker='*', color=sect_color,
    #     markersize=14, markeredgecolor='k')
    # ax.plot(dist[n_tn], ztop[n_tn] + 5, marker='o', color=sect_color,
    #     markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * fld_s
    # plot section
    cs = ax.pcolormesh(dist_se, zw_se, sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # labels
    ax.text(.6, .3, 'SECTION\nWA Coast', fontsize=fs, color='b',ha='center',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # # Tides
    # alpha = .4
    # ax = fig.add_subplot(436)
    # ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
    #     lw=3, alpha=alpha)
    # # time marker
    # ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
    #     marker='o', color='r', markersize=7)
    # # labels
    # ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
    #     alpha=alpha, fontsize=fs, horizontalalignment='right')
    # ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
    #     alpha=alpha, fontsize=fs, horizontalalignment='right')
    # # limits
    # ax.set_xlim(0,365)
    # ax.set_ylim(0,1.5)
    # ax.set_axis_off()

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
    x_e = x
    y_e = y
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = \
        pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

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
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.99,.97,'range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
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

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, ztop+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], ztop[n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], ztop[n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * fld_s
    # plot section
    cs = ax.pcolormesh(dist_se, zw_se, sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
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
