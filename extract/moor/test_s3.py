"""
Code to explore mooring-extraction from history files in s3 buckets.
"""

# imports
from time import time
tt00 = time()

import os
from lo_tools import Lfun
import xarray as xr
import s3fs

print('time for imports = %0.1f sec' % (time()-tt00))
# took 25 sec

gtx = 'cas7_t1_x11ab'
fstr = 'f2013.07.02'
fname = 'ocean_avg_0001.nc'

Ldir = Lfun.Lstart()

out_dir = Ldir['LOo'] / 'test_s3'
Lfun.make_dir(out_dir)
out_fn = out_dir / 'test.nc'

ENDPOINT = 'https://s3.kopah.uw.edu'
fs_s3 = s3fs.S3FileSystem(
    key=os.environ['AWS_ACCESS_KEY_ID'],
    secret=os.environ['AWS_SECRET_ACCESS_KEY'],
    client_kwargs={'endpoint_url': ENDPOINT},
)

s3_url = 's3://liveocean-pmacc/LO_roms/' + gtx + '/' + fstr + '/' + fname

tt0 = time()
s3_file_obj = fs_s3.open(s3_url, mode='rb',cache_type='blockcache', block_size=2**22)
print('time for fs_s3.open = %0.1f sec' % (time()-tt0))
# result: less than 1 sec

# this works
tt0 = time()
ds = xr.open_dataset(s3_file_obj, engine='h5netcdf')
print('time for xr.open_dataset = %0.1f sec' % (time()-tt0))
# result:  3-5 sec

tt0 = time()
a = ds.salt[0,:,10,10].values
print(a)
print('time to get salt profile = %0.1f sec' % (time()-tt0))
# result: 6 sec
