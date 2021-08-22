#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:55:45 2017

@author: PM5

Utility functions for pgrid.
"""
import numpy as np
import zfun
import h5py
import matfun
import netCDF4 as nc
import seawater as sw

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zrfun
import Lfun

def simple_grid(aa, res):
    dlat = aa[3] - aa[2]
    dlon = aa[1] - aa[0]
    mean_lat = np.mean(aa[2:])
    earth_rad = zfun.earth_rad(mean_lat)
    dlatm = earth_rad * np.pi * dlat/180
    dlonm = earth_rad * np.cos(np.pi*mean_lat/180) * np.pi * dlon/180
    nx = int(np.ceil(dlonm/res))
    ny = int(np.ceil(dlatm/res))
    plon_vec = np.linspace(aa[0], aa[1], nx)
    plat_vec = np.linspace(aa[2], aa[3], ny)
    return plon_vec, plat_vec

def stretched_grid(lon_list, x_res_list, lat_list, y_res_list):
    """
    Input:
    The four lists (units = degrees and meters) are points that define
    segmented vectors along which the resolution changes linearly.
    """
    plon_list = []
    plat_list = []
    if (len(lon_list) != len(x_res_list)) or (len(lat_list) != len(y_res_list)):
        print('Lists must be the same length')
        return np.array(plon_list), np.array(plat_list)

    lon_vec = np.array(lon_list)
    x_res_vec = np.array(x_res_list)
    lat_vec = np.array(lat_list)
    y_res_vec = np.array(y_res_list)

    R = zfun.earth_rad(np.mean(lat_vec))
    clat = np.cos(np.pi*np.mean(lat_vec)/180)

    lon = lon_list[0]
    plon_list.append(lon)
    while lon <= lon_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lon]), lon_vec)
        xres = (1-fr)*x_res_vec[i0] + fr*x_res_vec[i1]
        dlon = 180 * xres / (clat * R * np.pi)
        lon = lon + dlon
        plon_list.append(lon[0])
    lat = lat_list[0]
    plat_list.append(lat)
    while lat <= lat_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lat]), lat_vec)
        yres = (1-fr)*y_res_vec[i0] + fr*y_res_vec[i1]
        dlat = 180 * yres / (R * np.pi)
        lat = lat + dlat
        plat_list.append(lat[0])
    return np.array(plon_list), np.array(plat_list)

def load_bathy_nc(t_fn):
    ds = nc.Dataset(t_fn)
    tlon_vec = ds['lon'][:]
    tlat_vec = ds['lat'][:]
    tz = ds['z'][:]
    ds.close()
    return tlon_vec, tlat_vec, tz

def load_bathy(t_fn):
    try:
        a = h5py.File(t_fn)
        aa = dict()
        print(' using h5py')
        for item in a.keys():
            aa[item] = a[item][:]
            #print('    ' + item + ' ' + str(aa[item].shape))
        # reshaping because of how things are packed
        tlon_vec = a['lon'][:].flatten()
        tlat_vec = a['lat'][:].flatten()
        tz = a['z'][:].T
        #tmask = a['mask'][:].T
        # The field "mask" was created in the topo processing
        # as a matrix of ints with 1=ocean, 0=land.
        # We will supersede it here with our own masking.
        a.close()
    except:
        # for some reason h5py does not work with, for example,
        # psdem/PS_183m.mat
        tmat = matfun.loadmat(t_fn)
        print(' using matfun')
        # reshaping because of how things are packed
        tlon_vec = tmat['lon'][:].flatten()
        tlat_vec = tmat['lat'][:].flatten()
        tz = tmat['z'][:]
    return tlon_vec, tlat_vec, tz
    
def load_bathy2(t_fn, lon_vec, lat_vec):
    # load a section of the new NetCDF Smith-Sandwell bathy
    ds = nc.Dataset(t_fn)
    Lon = ds['lon'][:]
    Lat = ds['lat'][:]
    i0 = zfun.find_nearest_ind(Lon, lon_vec[0])
    i1 = zfun.find_nearest_ind(Lon, lon_vec[-1])
    j0 = zfun.find_nearest_ind(Lat, lat_vec[0])
    j1 = zfun.find_nearest_ind(Lat, lat_vec[-1])
    tlon_vec = Lon[i0-2:i1+3]
    tlat_vec = Lat[j0-2:j1+3]
    tz = ds['z'][j0-2:j1+3, i0-2:i1+3]
    ds.close()
    return tlon_vec, tlat_vec, tz

def make_nc(out_fn, plon, plat, lon, lat, z, m, dch):
    # get rid of old version
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    
    # create new NetCDF file
    foo = nc.Dataset(out_fn, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions
    M, L = lon.shape # use ROMS teminology
    size_dict = {'rho': (M, L),
                 'u': (M, L-1),
                 'v': (M-1, L),
                 'psi': (M-1, L-1),
                 'psi_ex': (M+1, L+1)}
    tag_list = ['rho', 'u', 'v', 'psi', 'psi_ex']
    for tag in tag_list:
        foo.createDimension('eta_'+tag, size_dict[tag][0])
        foo.createDimension('xi_'+tag, size_dict[tag][1])
    
    # create variables
    lon_var = dict()
    lat_var = dict()
    mask_var = dict()
    for tag in tag_list:
        lat_var[tag] = foo.createVariable('lat_'+tag, float, ('eta_'+tag, 'xi_'+tag))
        lon_var[tag] = foo.createVariable('lon_'+tag, float, ('eta_'+tag, 'xi_'+tag))
        mask_var[tag] = foo.createVariable('mask_'+tag, float, ('eta_'+tag, 'xi_'+tag))
    h_var = foo.createVariable('h', float, ('eta_rho', 'xi_rho'))
    f_var = foo.createVariable('f', float, ('eta_rho', 'xi_rho'))
    pm_var = foo.createVariable('pm', float, ('eta_rho', 'xi_rho'))
    pn_var = foo.createVariable('pn', float, ('eta_rho', 'xi_rho'))
    xl_var = foo.createVariable('xl', float)
    el_var = foo.createVariable('el', float)
    sph_var = foo.createVariable('spherical', 'c')
    
    # create other grids
    lon_dict = dict()
    lat_dict = dict()
    
    lon_dict['rho'] = lon
    lat_dict['rho'] = lat
    
    lon_dict['psi_ex'] = plon
    lat_dict['psi_ex'] = plat
    
    lon_dict['u'] = lon[:, :-1] + np.diff(lon, axis=1)/2
    lat_dict['u'] = lat[:, :-1]
    
    lon_dict['v'] = lon[:-1, :]
    lat_dict['v'] = lat[:-1, :] + np.diff(lat, axis=0)/2
    
    lon_dict['psi'] = plon[1:-1, 1:-1]
    lat_dict['psi'] = plat[1:-1, 1:-1]
    
    # create dx, dy and pm, pn
    ulon = plon[:-1, :]
    vlat = plat[:, :-1]
    R = zfun.earth_rad(np.mean(plat[:,1]))
    dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(ulon, axis=1)/180)
    dy = R * (np.pi*np.diff(vlat, axis=0)/180)
    
    # add data to fields
    for tag in tag_list:
        lon_var[tag][:] = lon_dict[tag]
        lat_var[tag][:] = lat_dict[tag]
        # start with all ones (unmasked for ROMS)
        mask_var[tag][:] = np.ones_like(lon_dict[tag], dtype=int)
    
    h_var[:] = -z
    pm_var[:] = 1/dx
    pn_var[:] = 1/dy
    f_var[:] = sw.f(lat_dict['rho'])
    xl_var[:] = dx[0,:].sum()
    el_var[:] = dy[:,0].sum()
    sph_var[:] = 'T'
    
    foo.close()
    
def make_nudgcoef(dch, out_dir):
    # Using info from From https://www.myroms.org/projects/src/ticket/627
   
    gfn = out_dir + 'grid.nc'
    sfn = out_dir + 'S_COORDINATE_INFO.csv'
    
    G = zrfun.get_basic_info(gfn, only_G=True)
    S_info_dict = Lfun.csv_to_dict(sfn)
    S = zrfun.get_S(S_info_dict)
    
    fn = out_dir + 'nudgcoef.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(fn, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions
    foo.createDimension('s_rho', S['N'])
    for tag in ['rho', 'u', 'v']:
        foo.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
        foo.createDimension('xi_'+tag, G['lon_'+tag].shape[1])
    
    fld2 = np.zeros_like(G['lon_rho'])
    nn = 6 # size of nudging band in gridpoints
    days_short = dch['nudging_days'][0]
    days_long = dch['nudging_days'][1]
    # make inverse time scales (d-1)
    t0 = 1/days_short
    t1 = 1/days_long
    
    if 'north' in dch['nudging_edges']:            
        for i in range(G['L']):
            for j in range(G['M'] - nn, G['M']):
                jj = j - G['M'] + nn + 1
                tnud = t1 + jj*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'south' in dch['nudging_edges']:            
        for i in range(G['L']):
            for j in range(nn):
                tnud = t0 - j*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'west' in dch['nudging_edges']:
        for i in range(nn):
            for j in range(G['M']):
                tnud = t0 - i*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])
    if 'east' in dch['nudging_edges']:
        for i in range(G['L'] - nn, G['L']):
            for j in range(G['M']):
                ii = i - G['L'] + nn + 1
                tnud = t1 + ii*(t0-t1)/nn
                fld2[j, i] = np.max([tnud, fld2[j, i]])

        
    fld3 = np.zeros((S['N'], G['M'], G['L']))
    for i in range(S['N']):
        fld3[i, :, :] = fld2.copy()
    
    # add 2D field data
    vv = foo.createVariable('M2_NudgeCoef', float, ('eta_rho', 'xi_rho'))
    vv.long_name = '2D momentum inverse nudging coefficients'
    vv.units = 'day-1'
    vv[:] = fld2
    
    # add 3D field data
    vn_dict = {'M3_NudgeCoef': '3D momentum inverse nudging coefficients',
               'tracer_NudgeCoef': 'generic tracer inverse nudging coefficients',
               'temp_NudgeCoef': 'temp inverse nudging coefficients',
               'salt_NudgeCoef': 'salt inverse nudging coefficients'}
    for vn in vn_dict.keys():
        vv = foo.createVariable(vn, float, ('s_rho', 'eta_rho', 'xi_rho'))
        vv.long_name = vn_dict[vn]
        vv.units = 'day-1'
        vv[:] = fld3
        
    foo.close()

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
    HH = HH + shift
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
    HH = HH - shift
    return HH
    
