"""
This makes the ocn forcing files for the updated ROMS, including the banas-fennel bio fields.

Testing:

run make_forcing_main.py -g cas7 -r forecast -d [today's date] -f ocn03 -test True

python make_forcing_main.py -g cas7 -r forecast -d [today's date] -f ocn03 -test True > test.log &

2024.09.12 This code is based on ocn01. The main difference is that it uses a new URL
because of a change at HYCOM. Also we omit the ncks-based Plan A and instead only use
a modified version of Plan B which gest single times. We also have to combine extractions
of single variables.

2024.11.19 jx: updated based on ocn02. It uses the same URL from HYCOM but applied different filters for the 1-hourly ssh (Godin) and 3-hourly u,v,t,s (hanning n=24). I also removed the condition "forecast" to enable generating backfill ocn from HYCOM URL directly

Performance:

To do:
- Ofun: use gsw instead of seawater for potential temp calculation
- Ofun: use xarray instaed of netCDF4
- main: hollow out arrays unless start_type = new
- throughout: stop using masked arrays

"""
import netCDF4 as nc
from pathlib import Path
import sys
from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np

import shutil, os
import pickle
import pandas as pd

from lo_tools import Lfun, zfun, zrfun, Ofun_nc
import Ofun
import Ofun_bio

# defaults
planB = False
planC = False
add_CTD = False
do_bio = True

# defaults related to testing
verbose = False
testing_ncks = False
testing_fmrc = False
testing_planC = False

if Ldir['testing']:
    verbose = True
    from importlib import reload
    reload(Ofun)
    reload(Ofun_bio)
else:
    pass
    
# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# *** automate when to set add_CTD to True ***
if this_dt == datetime(2012,10,7):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True
    
# this is where all the pre-processed files will go
h_out_dir = out_dir / 'Data'
# This already exists because it is created by the initialization, but we make sure it is clean
# so that testing, which does not make a clean version, is more consistent with operational use,
# which does.
Lfun.make_dir(h_out_dir, clean=True)

#if (Ldir['run_type'] == 'forecast') and (testing_planC == False):  #jx
if testing_planC == False:
    # this either gets new hycom files, or sets planB to True,
    # and planB may in turn set planC to True
    
    # form list of days to get, datetimes
    nd_f = np.ceil(Ldir['forecast_days'])  # jx: nd_f = 3
    dt0 = this_dt - timedelta(days=2)      # jx: previous 2 days
    dt1 = this_dt + timedelta(days=int(nd_f) + 2) # jx: later 5 days
#    dt_list_full = []
#    dtff = dt0
#    while dtff <= dt1:
#        dt_list_full.append(dtff)
#        dtff = dtff + timedelta(days=1)  # daily
#%    
    dt_list_full_hr = []
    dtff_hr = dt0
    dt2 = this_dt + timedelta(days=int(nd_f) + 2)
    while dtff_hr <= dt2:
        dt_list_full_hr.append(dtff_hr)
        dtff_hr = dtff_hr + timedelta(hours=1)  # hourly
#%       
    dt_list_full_3hr = []
    dtff_3hr = dt0
    dt3 = this_dt + timedelta(days=int(nd_f) + 2)
    while dtff_3hr <= dt3:
        dt_list_full_3hr.append(dtff_3hr)
        dtff_3hr = dtff_3hr + timedelta(hours=3)  # 3-hrly
#%%
    planB = True
    
    # Plan B: use fmrc one day at a time
    if planB == True:
        print('**** Using planB ****')
        result_dict['note'] = 'planB'
#%%
        # get the indices for extraction using ncks
        got_indices_ssh = False
        for ntries in range(10):
            if got_indices_ssh == False:
                # try again
                print('\nget_indices: ntries = ' + str(ntries))
                try:
                    #ind_dicts, got_indices = Ofun.get_indices(h_out_dir, dt_list_full, verbose=verbose) #jx: get e.g., ssh_tyx.nc
                    ind_dicts_ssh, got_indices_ssh = Ofun.get_indices_ssh(h_out_dir,   dt_list_full_hr,  verbose=verbose)
                    ind_dicts_u,   got_indices_u   = Ofun.get_indices_u(  h_out_dir,   dt_list_full_3hr, verbose=verbose)
                    ind_dicts_v,   got_indices_v   = Ofun.get_indices_v(  h_out_dir,   dt_list_full_3hr, verbose=verbose)
                    ind_dicts_s,   got_indices_s   = Ofun.get_indices_s(  h_out_dir,   dt_list_full_3hr, verbose=verbose)
                    ind_dicts_t,   got_indices_t   = Ofun.get_indices_t(  h_out_dir,   dt_list_full_3hr, verbose=verbose)
                except Exception as e:
                    print(e)
                    got_indices_ssh = False
            else:
                break
#%%            
#        if got_indices == True: # not use this
#            try:
#                for idt in range(len(dt_list_full)):
#                    got_fmrc = False # initialize each time
#                   data_out_fn =  h_out_dir / ('h' + dt_list_full[idt].strftime(Lfun.ds_fmt)+ '.nc')
#                    if verbose:
#                       print('\n' + str(data_out_fn))
#                    sys.stdout.flush()
#                    # get hycom forecast data from the web, and save it in the file "data_out_fn".
#                    # it tries 10 times before ending
#                    for ntries in range(10):
#                        if got_fmrc == False:
#                            # try again
#                            print('\nget_data_oneday: ntries = ' + str(ntries))
#                            got_fmrc = Ofun.get_data_oneday(idt, data_out_fn, ind_dicts, testing_fmrc, verbose=verbose)
#                        else:
#                            break
#                    if got_fmrc == False:
#                        # this should break out of the dtff loop at the first failure
#                        # and send the code to Plan C
#                        print('- error getting forecast files using fmrc')
#                        planC = True
#                        break
#            except Exception as e:
#                print(e)
#                planC = True
                
#%%------------  get ssh data ------------
        if got_indices_ssh == True:
            try:
                got_fmrc_ssh = False # initialize each time
                data_out_fn =  h_out_dir / ('h_ssh.nc') # get all hourly ssh at the same time
                if verbose:
                    print('\n' + str(data_out_fn))
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                for ntries in range(10):
                    if got_fmrc_ssh == False:
                        # try again
                        print('\nget_data_oneday: ntries = ' + str(ntries))
                        got_fmrc_ssh = Ofun.get_data_hrly(data_out_fn, ind_dicts_ssh, testing_fmrc, verbose=verbose)
                    else:
                        break
                if got_fmrc_ssh == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
            except Exception as e:
                print(e)
                planC = True
        elif got_indices_ssh == False:  # send the code to Plan C
            print('- error getting indices for forecast')
            planC = True
            
#%------------  get u data ------------
        if got_indices_u == True:
            try:
                got_fmrc_u = False # initialize each time
                data_out_fn =  h_out_dir / ('h_u.nc') # get all hourly ssh at the same time
                if verbose:
                    print('\n' + str(data_out_fn))
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                for ntries in range(10):
                    if got_fmrc_u == False:
                        # try again
                        print('\nget_data_oneday: ntries = ' + str(ntries))
                        got_fmrc_u = Ofun.get_data_hrly(data_out_fn, ind_dicts_u, testing_fmrc, verbose=verbose)
                    else:
                        break
                if got_fmrc_u == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
            except Exception as e:
                print(e)
                planC = True
        elif got_indices_u == False: # send the code to Plan C
            print('- error getting indices for forecast')
            planC = True

#%------------  get v data ------------
        if got_indices_v == True:
            try:
                got_fmrc_v = False # initialize each time
                data_out_fn =  h_out_dir / ('h_v.nc') # get all hourly ssh at the same time
                if verbose:
                    print('\n' + str(data_out_fn))
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                for ntries in range(10):
                    if got_fmrc_v == False:
                        # try again
                        print('\nget_data_oneday: ntries = ' + str(ntries))
                        got_fmrc_v = Ofun.get_data_hrly(data_out_fn, ind_dicts_v, testing_fmrc, verbose=verbose)
                    else:
                        break
                if got_fmrc_v == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
            except Exception as e:
                print(e)
                planC = True
        elif got_indices_v == False:
            # send the code to Plan C
            print('- error getting indices for forecast')
            planC = True

#%------------  get s data ------------
        if got_indices_s == True:
            try:
                got_fmrc_s = False # initialize each time
                data_out_fn =  h_out_dir / ('h_s.nc') # get all hourly ssh at the same time
                if verbose:
                    print('\n' + str(data_out_fn))
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                for ntries in range(10):
                    if got_fmrc_s == False:
                        # try again
                        print('\nget_data_oneday: ntries = ' + str(ntries))
                        got_fmrc_s = Ofun.get_data_hrly(data_out_fn, ind_dicts_s, testing_fmrc, verbose=verbose)
                    else:
                        break
                if got_fmrc_s == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
            except Exception as e:
                print(e)
                planC = True
        elif got_indices_s == False:
            # send the code to Plan C
            print('- error getting indices for forecast')
            planC = True

#%------------  get t data ------------
        if got_indices_t == True:
            try:
                got_fmrc_t = False # initialize each time
                data_out_fn =  h_out_dir / ('h_t.nc') # get all hourly ssh at the same time
                if verbose:
                    print('\n' + str(data_out_fn))
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                for ntries in range(10):
                    if got_fmrc_t == False:
                        # try again
                        print('\nget_data_oneday: ntries = ' + str(ntries))
                        got_fmrc_t = Ofun.get_data_hrly(data_out_fn, ind_dicts_t, testing_fmrc, verbose=verbose)
                    else:
                        break
                if got_fmrc_t == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
            except Exception as e:
                print(e)
                planC = True
        elif got_indices_t == False:
            # send the code to Plan C
            print('- error getting indices for forecast')
            planC = True
            
#%%
if planC == False:
    # process the hycom files, going from the original NetCDF extractions
    # to the processed pickled dicts.
    
#    if Ldir['run_type'] == 'backfill':
#        # Make a list of files to use from the hycom archive.
#        # This is a list of strings corresponding to NetCDF files
#        hnc_short_list = Ofun.get_hnc_short_list(this_dt, Ldir)
#        # step through those days and convert them to the same format
#        # of pickled dicts as used by the forecast
#        for fn in hnc_short_list:
#            a = Ofun.convert_extraction_oneday(fn)
#            dts = datetime.strftime(a['dt'], Lfun.ds_fmt)
#            out_fn = h_out_dir / ('h' + dts + '.p')
#            pickle.dump(a, open(out_fn, 'wb'))
#            
#    elif Ldir['run_type'] == 'forecast':
#        hnc_list = sorted([item.name for item in h_out_dir.iterdir()
#                if item.name[0]=='h' and item.name[-3:]=='.nc'])
#        #hnc_list.remove('h_ssh.nc')
#        for item in hnc_list:
#            a = Ofun.convert_extraction_oneday(h_out_dir / item) # get "h2024.09.14.p"
#            out_fn = h_out_dir / item.replace('.nc','.p')
#            pickle.dump(a, open(out_fn, 'wb'))
#            sys.stdout.flush()
            
#    hp_list = sorted([item.name for item in h_out_dir.iterdir()
#            if (item.name[0]=='h' and item.name[-2:]=='.p')])
        
#    # copy in the coordinates (assume those from first file work)
#    this_h_dict = pickle.load(open(h_out_dir / hp_list[0], 'rb'))
#    coord_dict = dict()
#    for vn in ['lon', 'lat', 'z']:
#        coord_dict[vn] = this_h_dict[vn]
#    pickle.dump(coord_dict, open(h_out_dir / 'coord_dict.p', 'wb'))
    
#    # filter in time
#    Ofun.time_filter(h_out_dir, hp_list, h_out_dir, Ldir) # get "fh2024.09.14.p"

#% ------------------------------------------------------
    # generate coor_dict.p
    ds = nc.Dataset(h_out_dir / 'h_s.nc')
    lon = ds.variables['lon'][:] - 360  # convert from 0:360 to -360:0 format
    lat = ds.variables['lat'][:]
    depth = ds.variables['depth'][:]
    z = -depth[::-1]
    coord_dict = dict()
    coord_dict['lon'] = lon
    coord_dict['lat'] = lat
    coord_dict['z'] = z
    pickle.dump(coord_dict, open(h_out_dir / 'coord_dict.p', 'wb'))

#% -------------- filter HYCOM fields -------------------
    # filter hourly ssh
    ds = xr.open_dataset(h_out_dir / 'h_ssh.nc') 
    time_ssh = ds.time.values
    ssh_lp = zfun.lowpass(ds.surf_el.values, f='godin')
    ds.close()
    # filter 3-hourly u
    ds = xr.open_dataset(h_out_dir / 'h_u.nc') 
    time_u = ds.time.values
    u3d_lp = zfun.lowpass(ds.water_u.values, f='hanning', n=24) # 3-day hanning
    u3d_lp = u3d_lp[:, ::-1, :,:]  # pack bottom to top
    ds.close()
    # filter 3-houry v
    ds = xr.open_dataset(h_out_dir / 'h_v.nc')
    time_v = ds.time.values
    v3d_lp = zfun.lowpass(ds.water_v.values, f='hanning', n=24)
    v3d_lp = v3d_lp[:, ::-1, :,:]  # pack bottom to top
    ds.close()
    # filter 3-hourly t
    ds = xr.open_dataset(h_out_dir / 'h_t.nc')
    time_t = ds.time.values
    t3d_lp = zfun.lowpass(ds.water_temp.values, f='hanning', n=24)
    t3d_lp = t3d_lp[:, ::-1, :,:]  # pack bottom to top
    ds.close()
    # filter 3-hourly s
    ds = xr.open_dataset(h_out_dir / 'h_s.nc')
    time_s = ds.time.values
    s3d_lp = zfun.lowpass(ds.salinity.values, f='hanning', n=24)
    s3d_lp = s3d_lp[:, ::-1, :,:]  # pack bottom to top
    ds.close()
    
#% ------------------------------------------------
    # save filtered data at daily interval
    # for example: fh2024.11.19.p, fh2024.11.20.p, fh2024.11.21.p, fh2024.11.22.p
    dts0 = Ldir['date_string']  # e.g. '2024.11.19'
    dt0 = datetime.strptime(dts0, '%Y.%m.%d')
    dt1 = dt0 + timedelta(days=3)
    dt = dt0
    while dt <= dt1:
       # print(dt)
        dts = datetime.strftime(dt, '%Y.%m.%d')
        out_name = 'fh' + dts + '.p'
        aa = dict()
        aa['dt'] = dt
        ix = np.argmin(np.abs(time_ssh - np.datetime64(dt)))
        aa['ssh'] = ssh_lp[ix,:,:]
        ix = np.argmin(np.abs(time_u - np.datetime64(dt)))
        aa['u3d'] = u3d_lp[ix,:,:,:]
        ix = np.argmin(np.abs(time_v - np.datetime64(dt)))
        aa['v3d'] = v3d_lp[ix,:,:,:]
        ix = np.argmin(np.abs(time_t - np.datetime64(dt)))
        aa['t3d'] = t3d_lp[ix,:,:,:]
        ix = np.argmin(np.abs(time_s - np.datetime64(dt)))
        aa['s3d'] = s3d_lp[ix,:,:,:]
    
        pickle.dump(aa, open(h_out_dir / out_name, 'wb'))
        dt += timedelta(days=1)
#%-------------------------------------------------------------
    
    # extrapolate
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(h_out_dir)
    fh_list = sorted([item.name for item in h_out_dir.iterdir()
            if item.name[:2]=='fh'])
    for fn in fh_list:
        print('-Extrapolating ' + fn)
        in_fn = h_out_dir / fn
        V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir)
        pickle.dump(V, open(h_out_dir / ('x' + fn), 'wb'))

    # and interpolate to ROMS format
    # get grid and S info
    G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
    S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
    S = zrfun.get_S(S_info_dict)
    # make list of files to process
    xfh_list = sorted([item.name for item in h_out_dir.iterdir()
            if item.name[:3]=='xfh'])
    # HyCOM grid info
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(h_out_dir)
    # load a dict of hycom fields
    dt_list = []
    count = 0
    c_dict = dict()
    zinds = Ofun.get_zinds(G['h'], S, z)
    for fn in xfh_list:
        print('-Interpolating ' + fn + ' to ROMS grid')
        b = pickle.load(open(h_out_dir / fn, 'rb'))
        dt_list.append(b['dt'])
        c = Ofun.get_interpolated(G, S, b, lon, lat, z, N, zinds)
        c_dict[count] = c
        count += 1
        
    # Write to ROMS forcing files
    # The fields we want to write are in c_dict, whose keys are indices of time.
    # Then c_dict[i] is another dict of arrays at that time, with keys:
    # ['ssh', 'ubar', 'vbar', 'theta', 's3d', 'u3d', 'v3d']
    # and the associated time vector is dt_list (datetimes).
    
    # Here I will write these to the dict V as expected below
    hycom_names = ['ssh', 'ubar', 'vbar', 'theta', 's3d', 'u3d', 'v3d']
    roms_names = ['zeta', 'ubar', 'vbar', 'temp', 'salt', 'u', 'v']
    names_dict = dict(zip(hycom_names, roms_names))
    
    # get sizes
    NZ = S['N']; NR = G['M']; NC = G['L']

    # Make the time vector.
    ot_vec = np.array([Lfun.datetime_to_modtime(item) for item in dt_list])
    NT = len(ot_vec)

    # Create a dict of fields for the state variables.
    V = dict()
    V['zeta'] = np.zeros((NT, NR, NC))
    V['ubar'] = np.zeros((NT, NR, NC-1))
    V['vbar'] = np.zeros((NT, NR-1, NC))
    V['salt'] = np.zeros((NT, NZ, NR, NC))
    V['temp'] = np.zeros((NT, NZ, NR, NC))
    V['u'] = np.zeros((NT, NZ, NR, NC-1))
    V['v'] = np.zeros((NT, NZ, NR-1, NC))
    
    # Fill the V dict
    for ii in range(NT):
        C = c_dict[ii]
        for vnh in hycom_names:
            vnr = names_dict[vnh]
            # note that the : here represents all axes after 0
            # and that it retains the correct shape
            V[vnr][ii, :] = C[vnh]
            
    if add_CTD:
        tt00 = time()
        z_rho = zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
        for vn in ['salt','temp']:
            fld = V[vn].copy()
            V[vn] = Ofun_bio.fill_polygons(fld, vn, G, z_rho, Ldir)
        print(' - add_CTD task for salt and temp: %0.2f sec' % (time()-tt00))
            
            
    # Create masks
    mr2 = np.ones((NT, NR, NC)) * G['mask_rho'].reshape((1, NR, NC))
    mr3 = np.ones((NT, NZ, NR, NC)) * G['mask_rho'].reshape((1, 1, NR, NC))
    mu2 = np.ones((NT, NR, NC-1)) * G['mask_u'].reshape((1, NR, NC-1))
    mu3 = np.ones((NT, NZ, NR, NC-1)) * G['mask_u'].reshape((1, 1, NR, NC-1))
    mv2 = np.ones((NT, NR-1, NC)) * G['mask_v'].reshape((1, NR-1, NC))
    mv3 = np.ones((NT, NZ, NR-1, NC)) * G['mask_v'].reshape((1, 1, NR-1, NC))

    # Apply masks
    V['zeta'][mr2==0] = np.nan
    V['ubar'][mu2==0] = np.nan
    V['vbar'][mv2==0] = np.nan
    V['salt'][mr3==0] = np.nan
    V['temp'][mr3==0] = np.nan
    V['u'][mu3==0] = np.nan
    V['v'][mv3==0] = np.nan
        
    # add bio variables if needed
    tt0 = time()
    if do_bio:
        bvn_list = ['NO3', 'NH4', 'chlorophyll', 'phytoplankton', 'zooplankton',
                'LdetritusN', 'SdetritusN', 'LdetritusC', 'SdetritusC',
                'TIC', 'alkalinity', 'oxygen']
        salt = V['salt'].copy()
        for bvn in bvn_list:
            V[bvn] = Ofun_bio.create_bio_var(salt, bvn)
        print('- Add bio variables: %0.2f sec' % (time()-tt0))
            
        if add_CTD:
            tt00 = time()
            bvn_list_short = ['NO3','TIC', 'alkalinity', 'oxygen']
            z_rho = zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
            for bvn in bvn_list_short:
                fld = V[bvn].copy()
                V[bvn] = Ofun_bio.fill_polygons(fld, bvn, G, z_rho, Ldir)
                V[bvn][mr3==0] = np.nan
            print(' - add_CTD task for bio variables: %0.2f sec' % (time()-tt00))
            
    # Write climatology file making use of zrfun.get_varinfo().
    #
    # NOTE: 2023.11.18 I should try filling the interior of the clm arrays with
    # nan's to make them smaller (except for start_type = new).
    tt0 = time()
    out_fn = out_dir / 'ocean_clm.nc'
    out_fn.unlink(missing_ok=True)
    ds = xr.Dataset()    
    for vn in V.keys():
        # tt00 = time()
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        # print(' -- time to get varinfo: %0.2f sec' % (time()-tt00))
        tname = vinfo['time_name']
        dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
        ds[vn] = (dims, V[vn])
        ds[vn].attrs['units'] = vinfo['units']
        ds[vn].attrs['long_name'] = vinfo['long_name']
        # time coordinate
        ds[tname] = ((tname,), ot_vec)
        ds[tname].attrs['units'] = Lfun.roms_time_units
    # and save to NetCDF
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()
    print('- Write clm file: %0.2f sec' % (time()-tt0))
    sys.stdout.flush()

elif planC == True:
    print('**** Using planC ****')
    result_dict['note'] = 'planC'
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format=Lfun.ds_fmt)
    clm_yesterday = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + ds_yesterday) / Ldir['frc'] / 'ocean_clm.nc'
    clm_today = out_dir / 'ocean_clm.nc'
    
    try:
        # new cleaner method: use open_dataset, update, and save to a new name
        ds = xr.open_dataset(clm_yesterday, decode_times=False)
        tname_list = [item for item in ds.coords if 'time' in item]
        for tname in tname_list:
            ot_vec = ds[tname].values
            ot_vec[-1] += 86400
            ds.update({tname: (('ocean_time',), ot_vec)})
            ds[tname].attrs['units'] = Lfun.roms_time_units
        ds.to_netcdf(clm_today)
        ds.close()
    except Exception as e:
        print(e)

if Ldir['start_type'] == 'new':
    # Write initial condition file if needed
    tt0 = time()
    in_fn = out_dir / 'ocean_clm.nc'
    out_fn = out_dir / 'ocean_ini.nc'
    out_fn.unlink(missing_ok=True)
    Ofun_nc.make_ini_file(in_fn, out_fn)
    print('- Write ini file: %0.2f sec' % (time()-tt0))
    sys.stdout.flush()

# Write boundary file
tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc.make_bry_file(in_fn, out_fn)
print('- Write bry file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    print(ds)
    ds.close()

# Check results
if Ldir['start_type'] == 'new':
    nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
else:
    nc_list = ['ocean_clm.nc', 'ocean_bry.nc']
    
if Ldir['testing']:
    # open datasets to have a peek manually
    dsc = xr.open_dataset(out_dir / 'ocean_clm.nc', decode_times=False)
    if Ldir['start_type'] == 'new':
        dsi = xr.open_dataset(out_dir / 'ocean_ini.nc', decode_times=False)
    dsb = xr.open_dataset(out_dir / 'ocean_bry.nc', decode_times=False)
        
result_dict['result'] = 'SUCCESS'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'FAIL'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
