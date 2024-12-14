"""
Code to test getting the hycom indices, and seeing if their ranges are
suitable for a forecast.

run test_get_indices.py -g cas7 -r forecast -d [today's date] -f ocn03

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
verbose = True
testing_do_not_get_indices = False # Set this to True to not get the hycom indices,
testing_do_not_get_data = False # Set this to True to not get the hycom data,
# e.g. if you already have it and want to speed up testing
testing_planB = False

if Ldir['testing']:
    verbose = True
    from importlib import reload
    reload(Ofun)
    reload(Ofun_bio)
    testing_planB = True
else:
    pass
    
# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# *** automate when to set add_CTD to True ***
# if this_dt == datetime(2012,10,7):
#     print('WARNING: adding CTD data to extrapolation!!')
#     add_CTD = True
    
# this is where all the pre-processed files will go
h_out_dir = out_dir / 'Data'
# This already exists because it is created by the initialization, but we make sure it is clean
# so that testing, which does not make a clean version, is more consistent with operational use,
# which does. We only clean it out if testing_do_not_get_data/indices == False.
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

# look at the contents of ind_dicts
for k in ind_dicts.keys():
    print('\n'+ k)
    a = ind_dicts[k]
    print('it0=%d, it1=%d, ix0=%d, ix1=%d, iy0=%d, iy1=%d' % (a['it0'], a['it1'], a['ix0'], a['ix1'], a['iy0'], a['iy1'], ))
    # and look at the times in each coordintate file
    out_fn = h_out_dir / (k + '_tyx.nc')
    ds = xr.open_dataset(out_fn)
    dt00 = pd.Timestamp(ds.time[a['it0']].values)
    dt11 = pd.Timestamp(ds.time[a['it1']].values)
    print('%s to %s' % (dt00.strftime('%Y.%m.%d %H:%M:%S'), dt11.strftime('%Y.%m.%d %H:%M:%S')))
    ds.close()