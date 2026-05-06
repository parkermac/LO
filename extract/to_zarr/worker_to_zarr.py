"""
Code to convert a single netcdf file to zarr.
"""

import xarray as xr
from time import time
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import argparse
import sys
from lo_tools import Lfun

parser = argparse.ArgumentParser()
parser.add_argument('-his_str', type=str)
parser.add_argument('-gtx', type=str)
parser.add_argument('-ds0', type=str)
args = parser.parse_args()

time_format = '%Y.%m.%d %H:%M:%S'
time_str = datetime.now().strftime(time_format)

tt0 = time()

f_string = 'f' + args.ds0

hh = ('000' + args.his_str)[-4:]
cmd_list = ['s5cmd','cp',
    's3://'+bucket_name+'/LO_roms/'+args.gtx+'/'+f_string+'/ocean_his_' + hh + '.nc']
    '/var/tmp/']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()

in_fn = '/var/tmp/ocean_his_' + hh + '.nc'
ds = xr.open_dataset(in_fn)

out_fn = '/var/tmp/ocean_his_' + hh + '.zarr'
ds.to_zarr(out_fn)

cmd_list = ['s5cmd','cp','--acl','public-read',
    out_fn, 's3://'+bucket_name+'/LO_roms/'+args.gtx+'_zarr/'+f_string+'/']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()

print('%s: time for to_zarr = %0.1f sec' % (time_str, time()-tt0))
sys.stdout.flush()


