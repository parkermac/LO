"""
Code to convert a single netcdf file to zarr.
"""

import xarray as xr
from time import time
from datetime import datetime, timedelta

import argparse
import sys
from lo_tools import Lfun

parser = argparse.ArgumentParser()
parser.add_argument('-his_str', type=str)
args = parser.parse_args()

time_format = '%Y.%m.%d %H:%M:%S'
time_str = datetime.now().strftime(time_format)

hh = ('000' + args.his_str)[-4:]

in_fn = '/var/tmp/ocean_his_' + hh + '.nc'
out_fn = '/var/tmp/ocean_his_' + hh + '.zarr'

tt0 = time()

ds = xr.open_dataset(in_fn)
ds.to_zarr(out_fn)

print('%s: time for to_zarr = %0.1f sec' % (time_str, time()-tt0))
sys.stdout.flush()
