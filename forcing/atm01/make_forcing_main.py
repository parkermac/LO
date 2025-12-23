"""
This is the main program for making the ATM forcing file.

Test on mac in ipython:

test a forecast (wgh2 runs faster than cas7):
run make_forcing_main.py -g wgh2 -r forecast -d 2019.07.04 -f atm01 -test True

Hour 19 of 2014.06.18 was bad for d3 (nans in the gridded files) so it is a useful test:
run make_forcing_main.py -g wgh2 -r backfill -d 2014.06.18 -f atm01 -test True
I used this on my mac to test new planB edits, either using or skipping (by renaming)
the bad d3 file.

test a forecast that will go to planB:
run make_forcing_main.py -g wgh2 -r forecast -d 2019.07.05 -f atm01 -test True
or
run make_forcing_main.py -g wgh2 -r forecast -d 2019.07.04 -f atm01 -test True -test_planB True

NEW 2025.07.30: For start_type = forecast this writes the forcing files
to separate day folders.

NOTE: atm_fun.py calls the old seawater routines in order to use the
"dist" method for getting grid angles for vector wind rotation.
Eventually we should just do this by hand.

NOTE: Only effect of -test True is to be verbose.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import os
import time
import shutil
import xarray as xr
import numpy as np
from scipy.spatial import cKDTree
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun, zfun, zrfun
import atm_fun as afun

verbose = False
if Ldir['testing']:
    verbose = True

# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir_data = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
    ('f' + Ldir['date_string']) / Ldir['frc'] / 'Data'
# We use this as a scratch directory for processsed files. Then either move the results to the
# main folder, or split them up into the forecast days.

# Set where are files located, and other situational choices.
do_d3 = True
do_d4 = True
wrf_dir = Ldir['data'] / 'wrf' # the default
if 'apogee' in Ldir['lo_env']:
    wrf_dir = Path('/dat1/parker/LO_data/wrf')
    
# Create list of hours
if Ldir['run_type'] == 'backfill':
    hr_vec = range(25)
elif Ldir['run_type'] == 'forecast':
    hr_max = Ldir['forecast_days'] * 24
    hr_vec = range(hr_max + 1)
    
# Create lists of input files.  These will be the complete lists
# regardless of whether or not the files exist, composed of the full
# path to each file (Path objects).
d_str = Ldir['date_string'].replace('.','')
in_dir = wrf_dir / (d_str + '00')
d2_list = []
d3_list = []
d4_list = []
for hr in hr_vec:
    hr_str = ('0' + str(hr))[-2:]
    d2_list.append(in_dir / ('wrfout.ocean_d2.' + d_str + '00.f' + hr_str + '.0000'))
    d3_list.append(in_dir / ('wrfout.ocean_d3.' + d_str + '00.f' + hr_str + '.0000'))
    d4_list.append(in_dir / ('wrfout.ocean_d4.' + d_str + '00.f' + hr_str + '.0000'))
    
# Create dict that relates a d2 filename to a time index (used when writing to NetCDF)
d2i_dict = {}
for i, v in enumerate(d2_list):
    d2i_dict[v] = i
        
# Check for existence of files. If any d2 are missing then exit.
planB = False
for fn in d2_list:
    if not fn.is_file():
        print('** Missing file: ' + str(fn))
        planB = True
        break

if Ldir['test_planB'] == True:
    planB = True
    
if planB == False:
    # For d3 and d4 just make sure we have the first one, so that we can get the grid
    for fn in [d3_list[0]]:
        if not fn.is_file():
            print('** Missing file: ' + str(fn))
            do_d3 = False
    for fn in [d4_list[0]]:
        if not fn.is_file():
            print('** Missing file: ' + str(fn))
            do_d4 = False

    # Create vector of time, in model format (seconds since whenever)
    dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
    mod_time_list = []
    for hr in hr_vec:
        dt = dt0 + timedelta(days=hr/24)
        mod_time = Lfun.datetime_to_modtime(dt)
        mod_time_list.append(mod_time)
    mod_time_vec = np.array(mod_time_list)
    NT = len(mod_time_vec)

    # Get ROMS model grid that we will interpolate to
    # get grid info, and some sizes
    G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
    NR = G['M']; NC = G['L']
    lon = G['lon_rho']
    lat = G['lat_rho']
    
    # Get WRF grids and grid size information (e.g. dx2_km = 12.5)
    # Note: lat, lon are only in the first file of the day (hour zero)
    lon2, lat2, dx2_km = afun.get_wrf_grid(d2_list[0])
    if do_d3:
        try:
            lon3, lat3, dx3_km = afun.get_wrf_grid(d3_list[0])
        except:
            do_d3 = False
    if do_d4:
        try:
            # sometimes there are empty files
            lon4, lat4, dx4_km = afun.get_wrf_grid(d4_list[0])
        except:
            do_d4 = False
        
    outvar_list = afun.outvar_list
        
    # Initialize output dict with zeros
    omat = np.zeros((NT, NR, NC))
    D = dict()
    for vn in outvar_list:
        D[vn] = omat.copy()

    # Find index to trim Eastern part of wrf fields
    lon_max = lon[0,-1] # easternmost edge of ROMS grid
    imax2 = zfun.find_nearest_ind(lon2[0,:], lon_max + .5)
    lon2 = lon2[:,:imax2]
    lat2 = lat2[:, :imax2]
    if do_d3:
        imax3 = zfun.find_nearest_ind(lon3[0,:], lon_max + .5)
        lon3 = lon3[:,:imax3]
        lat3 = lat3[:, :imax3]
    if do_d4:
        imax4 = zfun.find_nearest_ind(lon4[0,:], lon_max + .5)
        lon4 = lon4[:,:imax4]
        lat4 = lat4[:, :imax4]

    # Prepare coordinate arrays for interpolation
    XY = np.array((lon.flatten(), lat.flatten())).T # shape is (NR*NC, 2)
    XY2 = np.array((lon2.flatten(), lat2.flatten())).T
    if do_d3:
        XY3 = np.array((lon3.flatten(), lat3.flatten())).T
    if do_d4:
        XY4 = np.array((lon4.flatten(), lat4.flatten())).T

    # get nearest neighbor trees to use with wrf grids to interpolate
    # values from the wrf grids onto the ROMS grid
    IM2 = cKDTree(XY2).query(XY); IM2 = IM2[1]
    if do_d3:
        IM3 = cKDTree(XY3).query(XY); IM3 = IM3[1]
    if do_d4:
        IM4 = cKDTree(XY4).query(XY); IM4 = IM4[1]
    
    # Find coordinate rotation matrices to translate wrf velocity from
    # wrf grid directions to ROMS standard E+, N+
    ca2, sa2 = afun.get_angle(lon2, lat2)
    if do_d3:
        ca3, sa3 = afun.get_angle(lon3, lat3)
    if do_d4:
        ca4, sa4 = afun.get_angle(lon4, lat4)

    # Get the d3 and d4 masks. We appear to be avoiding a strip on the N.  Why?
    # We use these to overwrite just the parts where we have higher resolution wrf data.
    if do_d3:
        plon3_poly = np.concatenate((lon3[0,4:],lon3[:-5,-1],lon3[-5,4::-1],lon3[:-5:-1,4]))
        plat3_poly = np.concatenate((lat3[0,4:],lat3[:-5,-1],lat3[-5,4::-1],lat3[:-5:-1,4]))
        M3 = afun.get_indices_in_polygon(plon3_poly, plat3_poly, lon, lat)
    if do_d4:
        plon4_poly = np.concatenate((lon4[0,4:],lon4[:-5,-1],lon4[-5,4::-1],lon4[:-5:-1,4]))
        plat4_poly = np.concatenate((lat4[0,4:],lat4[:-5,-1],lat4[-5,4::-1],lat4[:-5:-1,4]))
        M4 = afun.get_indices_in_polygon(plon4_poly, plat4_poly, lon, lat)
    
    # MAIN TASK: loop over all hours

    dall_list = zip(d2_list, d3_list, d4_list)
    # Check out help(zip) to see how this works.  It creates an interable
    # that returns tuples made sequentially from entries of the things you zipped.
    # Note that this always works because we made our lists synthetically without regard
    # for if the files existed.
    for fn2, fn3, fn4 in dall_list:
        if verbose:
            print('Working on ' + str(fn2).split('/')[-1] + ' and etc.')
    
        # flags to allow processing more files
        do_this_d3 = True
        do_this_d4 = True
    
        # if we are missing a d3 or d4 file then we don't work on it
        if not fn3.is_file():
            print(' - missing ' + str(fn3))
            do_this_d3 = False
        if not fn4.is_file():
            print(' - missing ' + str(fn4))
            do_this_d4 = False

        try:
            # On 2022.06.07 one of the d2 files showed up with a time index of length 0
            # (it should have been 1).  This is designed to respond to that specific error.
            # I could wrap the whole job in a try-except but I like to keep it more specific.
            ov2_dict = afun.gather_and_process_fields(fn2, imax2, ca2, sa2, outvar_list)
        except Exception as e:
            print('Error in gather and process d2 fields - going to Plan B')
            print(e)
            planB = True
            break

        ovi2_dict = afun.interp_to_roms(ov2_dict, outvar_list, IM2, NR, NC)
    
        if do_this_d3:
            try:
                ov3_dict = afun.gather_and_process_fields(fn3, imax3, ca3, sa3, outvar_list)
                ovi3_dict = afun.interp_to_roms(ov3_dict, outvar_list, IM3, NR, NC)
            except:
                print(' - could not process ' + str(fn3))
                do_this_d3 = False
    
        if do_this_d4:
            try:
                ov4_dict = afun.gather_and_process_fields(fn4, imax4, ca4, sa4, outvar_list)
                ovi4_dict = afun.interp_to_roms(ov4_dict, outvar_list, IM4, NR, NC)
            except:
                print(' - could not process ' + str(fn4))
                do_this_d4 = False
    
        # combine the grids
        ovc_dict = dict()
        for ovn in outvar_list:
            v2 = ovi2_dict[ovn]
            v = v2.copy()
            if do_this_d3:
                v3 = ovi3_dict[ovn]
                v[M3] = v3[M3]
            if do_this_d4:
                v4 = ovi4_dict[ovn]
                v[M4] = v4[M4]
            if np.sum(np.isnan(v)) > 0:
                print('** WARNING Nans in combined output ' + ovn)
                print('Going to planB part 1')
                planB = True
                break # break from outvar_list loop
            ovc_dict[ovn] = v

        if planB == True:
            print('Going to planB part 2')
            break # break from dall_list loop
    
        # save to dict
        tt = d2i_dict[fn2]
        for vn in outvar_list:
            this_field = ovc_dict[vn]
            #this_field[G['mask_rho']==0] = np.nan
            D[vn][tt,:,:] = this_field
    
    # combine verything into a single dataset and save to Data
    for vn in outvar_list:
        out_fn = out_dir_data / (vn + '.nc')
        out_fn.unlink(missing_ok=True)
        ds = xr.Dataset()
        vinfo = zrfun.get_varinfo(vn)
        tname =  vinfo['time_name']
        dims = (tname,) + vinfo['space_dims_tup']
        # take time derivative or rain, converting [mm]
        # to [kg m-2 s-1]
        if vn == 'rain':
            rmat = D[vn].copy()
            rrmat = np.ones(rmat.shape)
            rrmat[:-1,:,:] = np.diff(rmat, axis=0) / 3600
            rrmat[-1,:,:] = rrmat[-2,:,:]
            D[vn] = rrmat

        # # HACK 2025.12.23 for negative Qair that showed up in a small patch
        # # and runined the forecast.
        # if vn == 'Qair':
        #     qmat = D[vn].copy()
        #     qmat[qmat<30] = 30
        #     D[vn] = qmat
        # # END HACK
        
        ds[vn] = (dims, D[vn])
        ds[vn].attrs['units'] = vinfo['units']
        ds[vn].attrs['long_name'] = vinfo['long_name']
        # time coordinate
        ds[tname] = ((tname,), mod_time_vec)
        ds[tname].attrs['units'] = Lfun.roms_time_units
        ds[tname].attrs['long_name'] = 'ocean time'
        # and save to NetCDF
        Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
        ds.to_netcdf(out_fn, encoding=Enc_dict)
        ds.close()

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Then reorganize the results:
    # - backfill: move from Data folder to main folder
    # - forecast: split using ncks into separate day folders
    out_fn_list = [] # start a complete list of all files created
    if Ldir['run_type'] == 'backfill':
        in_dir = out_dir_data
        out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
            ('f' + Ldir['date_string']) / Ldir['frc']
        for vn in outvar_list:
            in_fn = in_dir / (vn + '.nc')
            out_fn = out_dir / (vn + '.nc')
            out_fn.unlink(missing_ok=True)
            shutil.move(in_fn, out_fn)
            # since we move the files, we automatically clean up the Data folder
            out_fn_list.append(out_fn)
    elif Ldir['run_type'] == 'forecast':
        try:
            in_dir = out_dir_data
            ds00 = Ldir['date_string']
            dt00 = datetime.strptime(ds00, Lfun.ds_fmt)
            for day in range(Ldir['forecast_days']):
                dt0 = dt00 + timedelta(days=day)
                dt1 = dt00 + timedelta(days=day+1)
                ds0 = datetime.strftime(dt0, format=Lfun.ds_fmt)
                ds1 = datetime.strftime(dt1, format=Lfun.ds_fmt)
                ds0_for_ncks = datetime.strftime(dt0, format='%Y-%m-%d')
                ds1_for_ncks = datetime.strftime(dt1, format='%Y-%m-%d')
                out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
                    ('f' + ds0) / Ldir['frc']
                Lfun.make_dir(out_dir)
                for vn in outvar_list:
                    in_fn = in_dir / (vn + '.nc')
                    out_fn = out_dir / (vn + '.nc')
                    out_fn.unlink(missing_ok=True)
                    out_fn_list.append(out_fn)
                    ds = xr.open_dataset(in_fn)
                    # start creating the command list
                    cmd_list = ['ncks']
                    for vn1 in ds.data_vars:
                        vnda = ds[vn1] # get the DataArray for each variable
                        vnda_dim_tup = vnda.dims
                        for item in vnda_dim_tup:
                            if 'time' in item:
                                time_name = item
                                cmd_list += ['-d',time_name+','+ds0_for_ncks+','+ds1_for_ncks]
                    cmd_list.append(str(in_fn))
                    cmd_list.append(str(out_fn))
                    #print(cmd_list)
                    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
                    stdout, stderr = proc.communicate()
                    if len(stderr) > 0:
                        print('Error using ncks for %s' % (vn))
                        planB = True
                        break
                    ds.close()
        except Exception as e:
            print('Error using ncks to split forecast forcing')
            print(e)
            planB = True
    # Cleaning up
    if planB == False:
        Lfun.make_dir(out_dir_data, clean=True)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if planB == True:
    # We make this an if instead of elif because planB might have been set to True
    # in the if planB == True section above.
    #
    # Work plan
    # - backfill: copy yesterday to today and update time index
    # - forecast: copy day 3 of yesterday's forecast to day 3 of
    #       today's forecast and update time index
    result_dict['note'] = 'planB'
    print('**** Using planB ****')
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
    if Ldir['run_type'] == 'backfill':
        in_dt = dt_today - timedelta(days=1)
        out_dt = dt_today
    elif Ldir['run_type'] == 'forecast':
        nfd = Ldir['forecast_days']
        in_dt = dt_today + timedelta(days=nfd-2)
        out_dt = dt_today + timedelta(days=nfd-1)
    out_fn_list = []

    in_dstr = datetime.strftime(in_dt, format=Lfun.ds_fmt)
    out_dstr = datetime.strftime(out_dt, format=Lfun.ds_fmt)
    
    B_in_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
        ('f' + in_dstr) / Ldir['frc']
    B_out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
        ('f' + out_dstr) / Ldir['frc']
    Lfun.make_dir(B_out_dir, clean=True)
    # make the Info directory
    Lfun.make_dir(B_out_dir / 'Info')

    
    outvar_list = afun.outvar_list
    for ovn in outvar_list:
        B_in_fn = B_in_dir / (ovn + '.nc')
        B_out_fn = B_out_dir / (ovn + '.nc')
        out_fn_list.append(B_out_fn)
        # use open_dataset, update time, and save to a new name
        B_ds = xr.open_dataset(B_in_fn, decode_times=False)
        tname_list = [item for item in B_ds.coords if 'time' in item]
        for tname in tname_list:
            ot_vec = B_ds[tname].values
            ot_vec += 86400
            B_ds.update({tname: ((tname,), ot_vec)})
            B_ds[tname].attrs['units'] = Lfun.roms_time_units
        B_ds.to_netcdf(B_out_fn)
        B_ds.close()
        
# -------------------------------------------------------

# debugging
if verbose:
    for fn in out_fn_list:
        if fn.is_file():
            if 'Pair' in fn.name:
                print(str(fn))
                ds = xr.open_dataset(fn)
                print(ds)
                print('')

# test for success
result_dict['result'] = 'SUCCESS'
for fn in out_fn_list:
    if fn.is_file():
        pass
    else:
       result_dict['result'] = 'FAIL'
       
# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
