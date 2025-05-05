"""
This makes the ocn forcing files for ROMS, including the banas-fennel bio fields.

It uses GLORYS (glorys, hereafter) fields instead of HYCOM. Since the glorys fields
are daily averages is greatly simplifies the processing.

Testing:

run make_forcing_main.py -g cas7 -r backfill -d 2017.01.01 -f ocnG00 -test True


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
import pandas as pd

from lo_tools import Lfun, zfun, zrfun, Ofun_nc
from lo_tools import glorys_functions as gfun
import Ofun
import Ofun_bio

# defaults
planB = False
add_CTD = False
do_bio = True

# defaults related to testing
verbose = True
testing_do_not_get_data = False # Set this to True to not get the hycom data,
    # e.g. if you already have it and want to speed up testing
testing_planB = False

if Ldir['testing']:
    verbose = True
    from importlib import reload
    reload(Ofun)
    reload(Ofun_bio)
    # testing_planB = True
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
    # This either gets and processes glorys files, or sets planB to True.

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

    if Ldir['run_type'] == 'forecast':
        ## Need to get fields from glorys site
        
    elif Ldir['run_type'] == 'backfill':
        ## See if we already have the needed files
        ## otherwise go to the glorys site and get them

    got_data = True

    ## Test to see if we got the data

    if got_data:
        planB = False
    if got_data == False:
        print('- error getting data')
        planB = True
else:
    # do this if testing_planB == True
    planB = True

if planB == False:
    try:
        ## process the glorys fields to ROMS format
    except Exception as e:
        print(e)
        print('- error processing fields')
        planB = True
            

        roms_names = ['zeta', 'ubar', 'vbar', 'temp', 'salt', 'u', 'v']
        
                
        if add_CTD:
            tt00 = time()
            z_rho = zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
            for vn in ['salt','temp']:
                fld = V[vn].copy()
                V[vn] = Ofun_bio.fill_polygons(fld, vn, G, z_rho, Ldir)
            print(' - add_CTD task for salt and temp: %0.2f sec' % (time()-tt00))
                
            
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
        print(e)
        print('- error somwhere while creating clm file')
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
