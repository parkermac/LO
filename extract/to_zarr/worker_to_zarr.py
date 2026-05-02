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
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t2_x11b
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
# select time period and frequency
parser.add_argument('-d', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-his_str', type=str, default='01')

# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument to extract_argfun.intro(): ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

tt0 = time()

in_fn = '/var/tmp/h_' + Ldir['his_str'] + '.nc'
out_fn = '/var/tmp/h_' + Ldir['his_str'] + '.zarr'

ds = xr.open_dataset(in_fn)
ds.to_zarr(out_fn)

print('time for to_zarr = %0.1f sec' % (time()-tt0))
sys.stdout.flush()
