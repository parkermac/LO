#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:34:22 2017

@author: PM5

Module of functions for plotting grids in pgrid.
"""
import numpy as np
import pandas as pd

def get_plon_plat(using_old_grid, ds):
    if using_old_grid==True:
        # because older grids did not have lon,lat_psi_ex we create this
        # as an extension of lon,lat_psi
        plon0 = ds.variables['lon_psi'][:]
        plat0 = ds.variables['lat_psi'][:]
        dx = plon0[0,1] - plon0[0,0]
        dy = plat0[1,0] - plat0[0,0]
        ny0, nx0 = plon0.shape
        plon = np.nan * np.ones((ny0+2, nx0+2))
        plat = np.nan * np.ones((ny0+2, nx0+2))
        plon[1:-1, 1:-1] = plon0
        plat[1:-1, 1:-1] = plat0
        plon[:,0] = plon0[0,0] - dx
        plon[:,-1] = plon0[0,-1] + dx
        plon[0,:] = plon[1,:]
        plon[-1,:] = plon[-2,:]
        plat[0,:] = plat0[0,0] - dy
        plat[-1,:] = plat0[-1,0] + dy
        plat[:,0] = plat[:,1]
        plat[:,-1] = plat[:,-2]
    elif using_old_grid==False:
        plon = ds.variables['lon_psi_ex'][:]
        plat = ds.variables['lat_psi_ex'][:]
    return (plon, plat)
    
def get_grids(ds):
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
    return (lon_dict, lat_dict, mask_dict)
    
def add_river_tracks(Gr, ds, ax):
    # add river tracks and endpoints
    in_rfn = Gr['gdir'] + 'river_info.csv'
    try:
        df = pd.read_csv(in_rfn, index_col='rname')
    except FileNotFoundError:
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # get grids
    lon_dict, lat_dict, mask_dict = get_grids(ds)
    lonu = lon_dict['u']
    latu = lat_dict['u']
    lonv = lon_dict['v']
    latv = lat_dict['v']
    # plot river tracks
    for rn in df.index:
        fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        ax.plot(x, y, '-r', linewidth=2)
        ax.plot(x[-1], y[-1], '*r')    
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                    latu[row_dict_py[rn], col_dict_py[rn]], '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                    latu[row_dict_py[rn], col_dict_py[rn]], '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                    latv[row_dict_py[rn], col_dict_py[rn]], '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                    latv[row_dict_py[rn], col_dict_py[rn]], 'vb')
                    
def edit_mask_river_tracks(Gr, NR, ax):
    # add river tracks and endpoints for edit_mask.py
    in_rfn = Gr['gdir'] + 'river_info.csv'
    try:
        df = pd.read_csv(in_rfn, index_col='rname')
    except FileNotFoundError:
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # plot river tracks
    for rn in df.index:
        yy = NR - row_dict_py[rn] - 1
        xx = col_dict_py[rn]
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            ax.plot(xx+.5, yy, '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            ax.plot(xx+.5, yy, '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            ax.plot(xx, yy-.5, '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            ax.plot(xx, yy-.5, 'vb')    

def show_z_info(zm, ax):
    # find the max value of z (DEBUGGING)
    (rowmax, colmax) = np.unravel_index(np.argmax(zm), zm.shape)
    zmax = zm[rowmax, colmax]
    print('Max z = ' + str(zmax))
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    ax.plot(lon_rho[rowmax, colmax], lat_rho[rowmax, colmax], '*m', markersize=20)

def show_grids(ds, ax):
    lon_dict, lat_dict, mask_dict = get_grids(ds)
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))