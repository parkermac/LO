"""
This is the main program for making the OCN forcing file.  It uses HYCOM
extractions.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -f ocn0 -test True

This one adds CTD data for the initial condition:
run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2016.12.15 -f ocn0

Forecast version:
run make_forcing_main.py -g cas6 -t v3 -r forecast -s continuation -f ocn0 -d [TODAY]

NOTE: I haven't yet made much use of the Ldir['testing'] flag, and instead there are some ad hoc
flags like verbose, testing_ncks, an testing_fmrc.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import shutil, os
import pickle
import netCDF4 as nc
import numpy as np
import time
import pandas as pd

from lo_tools import Lfun, zrfun
import Ofun
import Ofun_nc
import Ofun_CTD
import Ofun_bio

if Ldir['testing']:
    from importlib import reload
    reload(Ofun)
    reload(Ofun_nc)
    reload(Ofun_CTD)
    reload(Ofun_bio)

# defaults
planB = False
planC = False
add_CTD = False
do_bio = True
verbose = False

if Ldir['testing']:
    verbose = True

testing_ncks = False
testing_fmrc = False

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# *** automate when to set add_CTD to True ***
if this_dt == datetime(2016,12,15):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True
    
# this is where all the pre-processed files will go
h_out_dir = out_dir / 'Data'
Lfun.make_dir(h_out_dir, clean=True)
    
if Ldir['run_type'] == 'forecast':
    # this either gets new hycom files, or sets planB to True,
    # and planB may in turn set planC to True
    
    # form list of days to get, datetimes
    nd_f = np.ceil(Ldir['forecast_days'])
    dt0 = this_dt - timedelta(days=2)
    dt1 = this_dt + timedelta(days=int(nd_f) + 2)
    dt_list_full = []
    dtff = dt0
    while dtff <= dt1:
        dt_list_full.append(dtff)
        dtff = dtff + timedelta(days=1)
    
    # Plan A: use ncks
    try:
        print('**** Using planA ****')
        result_dict['note'] = 'planA'
        got_ncks = False
        got_ncks = Ofun.get_data_ncks(h_out_dir, dt0, dt1, testing_ncks)
    except Exception as e:
        print(e)
        planB = True
    if got_ncks == False:
        print('- error getting forecast files using ncks')
        planB = True
    
    # Plan B: use fmrc one day at a time
    if planB == True:
        print('**** Using planB ****')
        result_dict['note'] = 'planB'
        try:
            for dtff in dt_list_full:
                got_fmrc = False
                data_out_fn =  h_out_dir / ('h' + dtff.strftime(Lfun.ds_fmt)+ '.nc')
                if verbose:
                    print('\n' + data_out_fn)
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_out_fn".
                # it tries 10 times before ending
                got_fmrc = Ofun.get_data_oneday(dtff, data_out_fn, testing_fmrc)
                if got_fmrc == False:
                    # this should break out of the dtff loop at the first failure
                    # and send the code to Plan C
                    print('- error getting forecast files using fmrc')
                    planC = True
                    break
        except Exception as e:
            print(e)
            planC = True
    
if planC == False:
    # process the hycom files, going from the original NetCDF extractions
    # to the processed pickled dicts.
    
    if Ldir['run_type'] == 'backfill':
        # Make a list of files to use from the hycom archive.
        # This is a list of strings corresponding to NetCDF files
        hnc_short_list = Ofun.get_hnc_short_list(this_dt, Ldir)
        # step through those days and convert them to the same format
        # of pickled dicts as used by the forecast
        for fn in hnc_short_list:
            a = Ofun.convert_extraction_oneday(fn)
            dts = datetime.strftime(a['dt'], Lfun.ds_fmt)
            out_fn = h_out_dir / ('h' + dts + '.p')
            pickle.dump(a, open(out_fn, 'wb'))
            
    elif Ldir['run_type'] == 'forecast':
        hnc_list = sorted([item.name for item in h_out_dir.iterdir()
                if item.name[0]=='h' and item.name[-3:]=='.nc'])
        for item in hnc_list:
            a = Ofun.convert_extraction_oneday(h_out_dir / item)
            out_fn = h_out_dir / item.replace('.nc','.p')
            pickle.dump(a, open(out_fn, 'wb'))
            sys.stdout.flush()
            
    hp_list = sorted([item.name for item in h_out_dir.iterdir()
            if (item.name[0]=='h' and item.name[-2:]=='.p')])
        
    # copy in the coordinates (assume those from first file work)
    this_h_dict = pickle.load(open(h_out_dir / hp_list[0], 'rb'))
    coord_dict = dict()
    for vn in ['lon', 'lat', 'z']:
        coord_dict[vn] = this_h_dict[vn]
    pickle.dump(coord_dict, open(h_out_dir / 'coord_dict.p', 'wb'))
    
    # filter in time
    Ofun.time_filter(h_out_dir, hp_list, h_out_dir, Ldir)

    # extrapolate
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(h_out_dir)
    fh_list = sorted([item.name for item in h_out_dir.iterdir()
            if item.name[:2]=='fh'])
    for fn in fh_list:
        print('-Extrapolating ' + fn)
        in_fn = h_out_dir / fn
        V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir,
            add_CTD=add_CTD)
        pickle.dump(V, open(h_out_dir / ('x' + fn), 'wb'))

    # and interpolate to ROMS format
    # get grid and S info
    G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
    S_fn = Ldir['grid'] / 'S_COORDINATE_INFO.csv'
    S_info_dict = pd.read_csv(S_fn, index_col='ITEMS').to_dict()['VALUES']
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
    Ofun_nc.make_clm_file(Ldir, out_dir, h_out_dir, c_dict, dt_list, S, G)
    
elif planC == True:
    print('**** Using planC ****')
    result_dict['note'] = 'planC'
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format=Lfun.ds_fmt)
    clm_yesterday = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + ds_yesterday) / Ldir['frc'] / 'ocean_clm.nc'
    clm_today = out_dir / 'ocean_clm.nc'
    shutil.copyfile(clm_yesterday, clm_today) # works with Path objects
    ds = nc.Dataset(clm_today, 'a')
    ot = ds['ocean_time'][:]
    ot[-1] += 86400
    for tname in ['ocean', 'salt', 'temp', 'v3d', 'v2d', 'zeta']:
        ds[tname + '_time'][:] = ot
    ds.close()

if do_bio and (planC==False):
    G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
    Ofun_bio.add_bio(out_dir, G, add_CTD=add_CTD)
    
Ofun_nc.make_ini_file(out_dir)
Ofun_nc.make_bry_file(out_dir)

# check results
nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
