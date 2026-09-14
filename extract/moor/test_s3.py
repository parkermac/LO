"""
Code to explore mooring-extraction from history files in s3 buckets.
"""

# imports
import sys, os
from lo_tools import Lfun, zrfun, zfun
import argparse
from time import time
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import numpy as np
import xarray as xr

import s3fs

gtx = 'cas7_t1_x11ab'
fstr = 'f2013.07.02'
fname = 'ocean_avg_0001.nc'

Ldir = Lfun.Lstart()

fn = Ldir['roms_out'] / gtx / fstr / fname

out_dir = Ldir['LOo'] / 'test_s3'
Lfun.make_dir(out_dir)
out_fn = out_dir / 'test.nc'

vn_list = 'h,zeta,salt,temp,u,v'

cmd_list = ['ncks',
    '-v', vn_list,
    '-d', 'xi_rho,'+str(10), '-d', 'eta_rho,'+str(10),
    '-d', 'xi_u,'+str(10), '-d', 'eta_u,'+str(10),
    '-d', 'xi_v,'+str(10), '-d', 'eta_v,'+str(10),
    '--mk_rec_dim','ocean_time']

# fs_s3 = s3fs.S3FileSystem(anon=False)

ENDPOINT = 'https://s3.kopah.uw.edu'
fs_s3 = s3fs.S3FileSystem(
    key=os.environ['AWS_ACCESS_KEY_ID'],
    secret=os.environ['AWS_SECRET_ACCESS_KEY'],
    client_kwargs={'endpoint_url': ENDPOINT},
)

#s3_url = 'https://s3.kopah.uw.edu/pm-share/kelly_2024.05.01_2024.05.31.nc'
#s3_url = 's3://pm-share/kelly_2024.05.01_2024.05.31.nc'
s3_url = 's3://liveocean-pmacc/LO_roms/' + gtx + '/' + fstr + '/' + fname

tt0 = time()
s3_file_obj = fs_s3.open(s3_url, mode='rb')
print('time for fs_s3.open = %0.1f sec' % (time()-tt0))

# this works
tt1 = time()
ds = xr.open_dataset(s3_file_obj, engine='h5netcdf',cache_type='blockcache', block_size=2**22)
print('time for xr.open_dataset = %0.1f sec' % (time()-tt1))

ds.salt[0,:,10,10].values

# tt_ncks = time()
# cmd_list += ['-O', str(fn), str(out_fn)]
# proc = Po(cmd_list, stdout=Pi, stderr=Pi)
# stdout, stderr = proc.communicate()
# print(' - time for ncks %0.2f sec' % (time()-tt_ncks))

# # check results
# ds = xr.open_dataset(out_fn)
# a = ds.salt.to_numpy().squeeze()
# print(a)