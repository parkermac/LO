"""
This is the main program for making a box extraction of the daily forecast and sending
the results to the APL server.

This is for for Susan Allen and Doug Latournell at UBC.

Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -r backfill -job ubc1 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's date string] -r forecast -job ubc1

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# +++ These are the main things to edit when changing jobs. +++

# Set the job name for LO/extract/box/extract_box.py, which refers to a choice
# in LO_user/extract/box/job_definitions.py.
box_job = 'ubc0'
# Note that box_job can be different from Ldir['job'] from the command line
# arguments.  Ldir['job'] controls the naming of where we archive the extraction:
# LO_output/post/[gtagex]/[fstring]/Ldir['job']

# This will be used for the file name that goes to the public server.
share_name = 'ubc'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# imports
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import shutil
from lo_tools import Lfun, zfun
import xarray as xr
import numpy as np
import pandas as pd

print(' - Creating extraction file for ' + Ldir['date_string'])

# create time range for box extraction
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (box_job + '_' + ds0 + '_' + ds1 + '.nc')

# this it the name of the file we will copy/move the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn_raw = out_dir / (share_name + '_raw.nc')
out_fn = out_dir / (share_name + '.nc')

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', ds0, '-1', ds1, '-lt', 'hourly', '-job', box_job]
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# copy/move the extraction to the expected "post" place
if Ldir['testing']:
    # when testing we keep the original extraction to make it easier to plot
    shutil.copyfile(out_fn0, out_fn_raw)
else:
    shutil.move(out_fn0, out_fn_raw)
print('\nPath to raw file:\n%s' % (str(out_fn_raw)))

# ===========================================================
# more processing to low-pass

dsr = xr.open_dataset(out_fn_raw)
ds = xr.Dataset()
# first add a weighting filter to dsr
NT = len(dsr.ocean_time)
if NT == 73:
    filt = zfun.godin_shape()
else:
    filt = zfun.hanning_shape(NT - 2)
filt = np.concatenate((np.array([0]),filt,np.array([0])))
dsr['filt'] = (('ocean_time'), filt)
# then make weighted means and write to ds
tt0 = time()
for vn in dsr.data_vars:
    ds[vn] = dsr[vn].weighted(dsr.filt).mean(dim='ocean_time', keep_attrs=True)
# The mean operation above removes the ocean_time dimension, which Susan would like,
# so we add it back in and then add the time as a coordinate.
ds = ds.expand_dims('ocean_time')
a = dsr.ocean_time.mean().values
ds['ocean_time'] = (('ocean_time'), pd.DatetimeIndex([a]))
# the coordinate step in the line above wasted about a day of my time until I finally tried
# making it a pandas DatetimeIndex object!!

print('Time to make weighted mean: %0.2f sec'% (time()- tt0))
sys.stdout.flush()
ds.to_netcdf(out_fn)

dsr.close()
ds.close()

# compress the resulting file
tt0 = time()
ds = xr.load_dataset(out_fn)
# need to load in order to do the compression
enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
Enc_dict = {vn:enc_dict for vn in ds.data_vars if 's_rho' in ds[vn].dims}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
print('Time to compress = %0.2f sec' % (time()- tt0))
# ===========================================================

# copy the file to the server
post_argfun.copy_to_server(Ldir, out_fn)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
