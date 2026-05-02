"""
Code to convert a single netcdf file to zarr.

Save to kopah?
"""

import xarray as xr
from time import time

import argparse
import sys
from lo_tools import Lfun

parser = argparse.ArgumentParser()
parser.add_argument('-his_str', type=str, default='01')
args = parser.parse_args()

tt0 = time()

in_fn = '/var/tmp/h_' + args.his_str + '.nc'
out_fn = '/var/tmp/h_' + args.his_str + '.zarr'

ds = xr.open_dataset(in_fn)
ds.to_zarr(out_fn)

print('time for to_zarr = %0.1f sec' % (time()-tt0))
sys.stdout.flush()
