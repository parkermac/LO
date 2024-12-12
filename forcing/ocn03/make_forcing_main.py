"""
This makes the ocn forcing files for the updated ROMS, including the banas-fennel bio fields.

Testing:

run make_forcing_main.py -g cas7 -r forecast -d [today's date] -f ocn03 -test True

python make_forcing_main.py -g cas7 -r forecast -d [today's date] -f ocn03 -test True > test.log &

2024.11.19 jx and PM: updated based on ocn02. The main improvement is that it handles the
fact that the new version of hycom (as of August 2024) has tides. This caused serious problems
with our old processing scheme which was based on daily snapshots. These aliased in a huge tidal
signal with fortnightly variation of ssh and velocity. This new version removes the tides by getting 1-hourly ssh
(and using a Godin filter) and 3-hourly u,v,t,s (and using a hanning n=24 filter - so three days).

Performance:
Because we download so much data this is slower, but still only takes about 15 minutes [CHECK]
for a 3-day forecast.

To do:
- main: hollow out arrays unless start_type = new
- throughout: stop using masked arrays

"""

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
add_CTD = False
do_bio = True

# defaults related to testing
verbose = False
testing_do_not_get_indices = False # Set this to True to not get the hycom indices,
testing_do_not_get_data = False # Set this to True to not get the hycom data,
# e.g. if you already have it and want to speed up testing
testing_planB = False

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
# which does. We only clean it out if testing_do_not_get_data == False.
if (testing_do_not_get_data == False) or (testing_do_not_get_indices == False):
    Lfun.make_dir(h_out_dir, clean=True)
else:
    print('WARNING: skipped extracting hycom data or getting indices')

if testing_planB == False:
    # This either gets and processes new hycom files, or sets planB to True.

    # Get the time limits for the hycom extraction. These go from the start of the
    # day to the end of the day (or day 3 for the forecast) with 2 days of padding
    # on either end to allow for filtering.
    if Ldir['run_type'] == 'forecast':
        nd_f = np.ceil(Ldir['forecast_days'])
        dt0 = this_dt - timedelta(days=2)      
        dt1 = this_dt + timedelta(days=int(nd_f) + 2)
    elif Ldir['run_type'] == 'backfill':
        dt0 = this_dt - timedelta(days=2)      
        dt1 = this_dt + timedelta(days=3)
    if verbose:
        print('dt0 = ' + str(dt0))
        print('dt1 = ' + str(dt1))

    # create lists and dicts of variable names and url's to hycom fields
    if Ldir['run_type'] == 'forecast':
        # New ones to try for forecast
        url_ssh  = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ssh/FMRC_ESPC-D-V02_ssh_best.ncd'
        url_uvel = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_u3z/FMRC_ESPC-D-V02_u3z_best.ncd'
        url_vvel = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_v3z/FMRC_ESPC-D-V02_v3z_best.ncd'
        url_temp = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_t3z/FMRC_ESPC-D-V02_t3z_best.ncd'
        url_salt = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_s3z/FMRC_ESPC-D-V02_s3z_best.ncd'
        # lists and dicts
    elif Ldir['run_type'] == 'backfill':
        # New ones to try for backfill
        url_ssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/ssh/2024'
        # url_Sssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/Sssh/2024'
        url_uvel  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/u3z/2024'
        url_vvel  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/v3z/2024'
        url_temp  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/t3z/2024'
        url_salt  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/s3z/2024'
        # lists and dicts
    hkeys = ['ssh','u', 'v','t', 's']
    url_list = [url_ssh, url_uvel, url_vvel, url_temp, url_salt]
    hycom_var_list = ["surf_el", "water_u", "water_v", "water_temp", "salinity"]
    url_dict = dict(zip(hkeys,url_list))
    hycom_var_dict = dict(zip(hkeys,hycom_var_list))

    # get the indices for extraction using ncks
    got_indices = False
    if testing_do_not_get_indices:
        got_indices = True
    else:
        for ntries in range(10):
            if got_indices == False:
                # try again
                print('\nget_indices: ntries = ' + str(ntries))
                try:
                    ind_dicts, got_indices = Ofun.get_indices(h_out_dir, dt0, dt1, url_dict, verbose=verbose)
                except Exception as e:
                    print(e)
                    got_indices = False
            else:
                break
                
    if got_indices == True:
        if testing_do_not_get_data:
            got_hycom_data = True
        else:
            try:
                for hkey in hkeys:
                    got_hycom_data = False # initialize each time
                    data_out_fn =  h_out_dir / ('h_' + hkey + '.nc') # get all hourly ssh at the same time
                    if verbose:
                        print('\n' + str(data_out_fn))
                    sys.stdout.flush()
                    # get hycom forecast data from the web, and save it in the file "data_out_fn".
                    # it tries 10 times before ending
                    for ntries in range(10):
                        if got_hycom_data == False:
                            # try again
                            print('\nget_hycom_data: ntries = ' + str(ntries))
                            got_hycom_data = Ofun.get_hycom_data(data_out_fn, hkey, ind_dicts, url_dict, hycom_var_dict, verbose=verbose)
                        else:
                            break
                    if got_hycom_data == False:
                        # this should break out of the dtff loop at the first failure
                        # and send the code to Plan B
                        print('- error getting hycom files')
                        planB = True
            except Exception as e:
                print(e)
                planB = True
    elif got_indices == False:  # send the code to Plan B
        print('- error getting indices')
        planB = True
else:
    # do this if testing_planB == True
    planB = True

if planB == False:
    try:
        # process the hycom files, going from the original NetCDF extractions
        # to the processed pickled dicts.
        
        # ------------------------------------------------------
        # generate coor_dict.p
        ds = xr.open_dataset(h_out_dir / 'h_s.nc')
        lon = ds.lon.values - 360  # convert from 0:360 to -360:0 format
        lat = ds.lat.values
        depth = ds.depth.values
        z = -depth[::-1] # repack bottom to top
        coord_dict = dict()
        coord_dict['lon'] = lon
        coord_dict['lat'] = lat
        coord_dict['z'] = z
        pickle.dump(coord_dict, open(h_out_dir / 'coord_dict.p', 'wb'))

        # -------------- filter HYCOM fields -------------------
        lp_time_dict = dict()
        lp_dict = dict()
        for hkey in hkeys:
            if hkey == 'ssh':
                # filter hourly ssh
                ds = xr.open_dataset(h_out_dir / ('h_' + hkey + '.nc')) 
                lp_time_dict[hkey] = ds.time.values
                lp_dict[hkey] = zfun.lowpass(ds[hycom_var_dict[hkey]].values, f='godin')
                ds.close()
            else:
                # filter 3-hourly 3D fields
                ds = xr.open_dataset(h_out_dir / ('h_' + hkey + '.nc')) 
                lp_time_dict[hkey] = ds.time.values
                this_lp = zfun.lowpass(ds[hycom_var_dict[hkey]].values, f='hanning', n=24) # 3-day hanning
                this_lp = this_lp[:, ::-1, :,:]  # pack bottom to top
                lp_dict[hkey] = this_lp
                ds.close()
        
        # ------------------------------------------------
        # save filtered data at daily interval
        # for example: fh2024.11.19.p, fh2024.11.20.p, fh2024.11.21.p, fh2024.11.22.p
        dt00 = this_dt
        if Ldir['run_type'] == 'forecast':
            nd_f = np.ceil(Ldir['forecast_days'])
            dt11 = dt00 + timedelta(days=int(nd_f))
        elif Ldir['run_type'] == 'backfill':
            dt11 = dt00 + timedelta(days=1)
        # We use a different naming convention at this point to be compatible with earlier processing steps
        # that we would like to re-use from ocn02.
        alt_name_dict = {'ssh':'ssh', 'u':'u3d', 'v':'v3d', 't':'t3d', 's':'s3d'}
        dt = dt00
        while dt <= dt11:
        # print(dt)
            dts = datetime.strftime(dt, Lfun.ds_fmt)
            out_name = 'fh' + dts + '.p'
            aa = dict()
            aa['dt'] = dt
            for hkey in hkeys:
                ix = np.argmin(np.abs(lp_time_dict[hkey] - np.datetime64(dt)))
                aa[alt_name_dict[hkey]] = lp_dict[hkey][ix,:] # the last : expands as needed
                # if verbose:
                #     print('checking on subsampling for %s %s' % (hkey,dts))
                #     print(aa[alt_name_dict[hkey]].shape)
            pickle.dump(aa, open(h_out_dir / out_name, 'wb'))
            dt += timedelta(days=1)

        #----------- back to original processing steps -------------------------
        
        # extrapolate (also get ubar and vbar and convert t to ptmp)
        lon, lat, z, L, M, N, X, Y = Ofun.get_coords(h_out_dir)
        fh_list = sorted([item.name for item in h_out_dir.iterdir()
                if item.name[:2]=='fh'])
        for fn in fh_list:
            print('-Extrapolating ' + fn)
            in_fn = h_out_dir / fn
            Vex = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir)
            # # debugging
            # for k in Vex.keys():
            #     if k == 'dt':
            #         pass
            #     else:
            #         print('%s %d' % (k, np.isnan(Vex[k]).sum()))
            #         print('Vex: max and min of ' + k)
            #         print(np.min(Vex[k]))
            #         print(np.max(Vex[k]))
            pickle.dump(Vex, open(h_out_dir / ('x' + fn), 'wb'))

        # # debugging
        # sys.exit()

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

    except Exception as e:
        # This try-except structure really should be more granular, with multiple
        # ones around each of the staps above.
        print(e)
        planB = True

elif planB == True:

    # planB means that we use the ocean_clm.nc file from the day before and change its last time value
    # to be a day later.
    print('**** Using planB ****')
    result_dict['note'] = 'planB'
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format=Lfun.ds_fmt)
    clm_yesterday = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + ds_yesterday) / Ldir['frc'] / 'ocean_clm.nc'
    clm_today = out_dir / 'ocean_clm.nc'
    try:
        # use open_dataset, update, and save to a new name
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
