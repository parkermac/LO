"""
This is a python job to do part of a tidal average.
"""
# imports
import sys
import argparse
from lo_tools import Lfun, zfun, zrfun
import os
from time import time
import numpy as np
import xarray as xr
from datetime import datetime, timedelta

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v0_live
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# NOTE: ds0 is day 1 of the three days required to make the tidally averaged file.
# The eventual product will have a timestamp at noon of day 2.
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-ii0', type=int) # index to start
parser.add_argument('-ii1', type=int) # index to end
parser.add_argument('-fnum', type=int) # number to use for filename
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
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
    
# output location
temp_out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'lowpass' / ('temp_' + Ldir['ds0'])
# Lfun.make_dir(temp_out_dir) # Assume it was created by the calling function.

# set up to do the partial lowpass
ds0 = Ldir['ds0']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=2)
ds1 = dt1.strftime(Lfun.ds_fmt)
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1) # length = 73
fn_list = fn_list[1:-1] # trim to length = 71

# filter shape
gs = zfun.godin_shape() # length 71, sum = 1

if Ldir['testing']:
    vn_list = ['h','zeta','salt']
    #vn_list = ['zeta','salt','temp','u','v']
else:
    vn_list = ['h','zeta','salt','temp','u','v','w',
    'NO3','phytoplankton','zooplankton',
    'detritus','Ldetritus','oxygen',
    'TIC','alkalinity',
    'Pair','Uwind','Vwind','shflux','ssflux','latent','sensible','lwrad','swrad',
    'sustr','svstr','bustr','bvstr']

for ii in range(Ldir['ii0'], Ldir['ii1']+1):
    ds = xr.open_dataset(fn_list[ii])
    if ii == Ldir['ii0']:
        lp = (ds[vn_list] * gs[ii]).squeeze(dim='ocean_time', drop=True).compute()
    else:
        lp = (lp + (ds[vn_list] * gs[ii]).squeeze(dim='ocean_time', drop=True)).compute()

# save the temp file to NetCDF
fstr = ('000' + str(Ldir['fnum']))[-2:]
lp.to_netcdf(temp_out_dir / ('lp_temp_' + fstr + '.nc'))

