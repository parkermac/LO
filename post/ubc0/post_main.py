"""
This is the main program for making the extraction for for Susan Allen
and Doug Latournell at UBC.

Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job ubc0

Test on apogee
run post_main.py -gtx cas6_v0_u0kb -ro 0 -d 2021.10.28 -job ubc0

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# imports
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import shutil
import xarray as xr
from lo_tools import zfun
import numpy as np

start_time = datetime.now()

print(' - Creating ubc file for ' + Ldir['date_string'])

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (Ldir['job'] + '_' + Ldir['date_string'] + '_' + Ldir['date_string'] + '.nc')

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn_raw = out_dir / 'UBC.nc'
out_fn = out_dir / 'low_pass_UBC.nc'
out_fn_raw.unlink(missing_ok=True)
out_fn.unlink(missing_ok=True)

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', Ldir['date_string'], '-1', Ldir['date_string'],
    '-lt', 'allhours', '-job', 'ubc0']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Time for initial extraction = %0.2f sec' % (time()-tt0))

# move the extraction to the expected "post" place
if Ldir['testing']:
    # if testing we keep the extraction around to look at
    shutil.copyfile(out_fn0, out_fn_raw)
else:
    shutil.move(out_fn0, out_fn_raw)

# and do the low pass
dsr = xr.load_dataset(out_fn_raw)
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
print('Time to make weighted mean: %0.2f sec'% (time()- tt0))
sys.stdout.flush()
ds.to_netcdf(out_fn)
dsr.close()
ds.close()

# squeeze and compress the resulting file
tt0 = time()
ds = xr.load_dataset(out_fn)
ds = ds.squeeze() # remove singleton dimensions
enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
Enc_dict = {vn:enc_dict for vn in ds.data_vars if 's_rho' in ds[vn].dims}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
print('Time to compress = %0.2f sec' % (time()- tt0))

# clean up
out_fn_raw.unlink(missing_ok=True)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
