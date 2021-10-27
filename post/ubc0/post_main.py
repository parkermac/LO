"""
This is the main program for Susan Allen and Doug Latournell at UBC.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job ubc0

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# recoded for LO, 2021.10.27 PM.

from lo_tools import Lfun, zfun, zrfun

from importlib import reload
import UBC_subdomain
reload(UBC_subdomain)

import os
import numpy as np
import shutil

# generate a list of all the files
in_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string'])
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir / ff) for ff in fn_list_raw if 'ocean_his' in ff]
fn_list.sort()
if len(fn_list) != 73:
    print('Warning: We have %d history files' % (len(fn_list)))
    
G = zrfun.get_basic_info(fn_list[0], only_G=True)

lon_vec = G['lon_rho'][0,:]
lat_vec = G['lat_rho'][:,0]

# Bounds of subdomain (rho grid) from Susan Allen 2018_11
lon0 = -125.016452048434
lon1 = -124.494612925929
lat0 = 48.3169463809796
lat1 = 48.7515055163539

# find indices containing these bonds
i0, i1, ifr = zfun.get_interpolant(np.array([lon0,lon1]), lon_vec)
j0, j1, jfr = zfun.get_interpolant(np.array([lat0,lat1]), lat_vec)
if np.nan in ifr:
    print('Warning: nan in ifr')
if np.nan in jfr:
    print('Warning: nan in jfr')
XBS = [i0[0], i1[1]]
YBS = [j0[0], j1[1]]
print(XBS)
print(YBS)

# the output file name
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ubc.nc'
out_fn.unlink(missing_ok=True)
    
# make a temporary place for the output
temp_dir = out_dir / 'ubc_temp'
Lfun.make_dir(temp_dir, clean=True)

# do the extraction for the whole list
UBC_subdomain.get_UBC_subdomain(fn_list, temp_dir, XBS, YBS)

# prepare for the low_pass
out_fn_list_raw = os.listdir(temp_dir)
out_fn_list = [(temp_dir / ff) for ff in out_fn_list_raw if 'ocean_his' in ff]
out_fn_list.sort()
# create the filter
out_fn_list_short = out_fn_list[1:-1] # trim to get to 71 for a forecast
nf = len(out_fn_list_short)
if nf == 71:
    print(' - Using Godin filter')
    filt0 = zfun.godin_shape()
else:
    print(' - Using Hanning filter for list length = ' + str(nf))
    filt0 = zfun.hanning_shape(nf)
# and do the low pass
zrfun.roms_low_pass(out_fn_list_short, out_fn, filt0, exclude=[])

# get rid of the temp directory
shutil.rmtree(temp_dir, ignore_errors=True)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
