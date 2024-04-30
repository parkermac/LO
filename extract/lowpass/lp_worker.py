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
# NOTE: dslp is day 2 of the three days required to make the tidally averaged file.
# The eventual product will have a timestamp at noon of dslp.
parser.add_argument('-dslp', type=str) # e.g. 2019.07.04
parser.add_argument('-ii0', type=int) # index to start
parser.add_argument('-ii1', type=int) # index to end
parser.add_argument('-fnum', type=int) # number to use for filename
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
parser.add_argument('-v', '--verbose', default=False, type=Lfun.boolean_string)
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
temp_out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'lowpass' / ('temp_' + Ldir['dslp'])
# Lfun.make_dir(temp_out_dir) # Assume it was created by the calling function.

# set up to do the partial lowpass
dslp = Ldir['dslp']
dt0 = datetime.strptime(dslp, Lfun.ds_fmt) - timedelta(days=1)
dt1 = dt0 + timedelta(days=2)
ds0 = dt0.strftime(Lfun.ds_fmt)
ds1 = dt1.strftime(Lfun.ds_fmt)
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1) # length = 73
fn_list = fn_list[1:-1] # trim to length = 71

# filter shape
gs = zfun.godin_shape() # length 71, sum = 1

# determine variables to process
ds = xr.open_dataset(fn_list[0])
vnl0 = list(ds.data_vars)
vnl1 = [vn for vn in vnl0 if 'time' in ds[vn].attrs.keys()]
ds.close()

# starter list
vn_list = ['zeta','salt','temp','u','v','w','ubar','vbar','AKs','AKv']

# add bio variables if available
if ('NO3' in vnl1) and ('NH4' not in vnl1): # old bio version
    vn_list += ['NO3','phytoplankton','zooplankton',
        'detritus','Ldetritus','oxygen','TIC','alkalinity']
elif ('NO3' in vnl1) and ('NH4' in vnl1): # new bio version
    vn_list += ['NO3','NH4','phytoplankton','zooplankton',
        'LdetritusN','SdetritusN','oxygen','TIC','alkalinity']

# add atm variables if available
if 'Pair' in vnl1:
    vn_list += ['Pair','Uwind','Vwind','shflux','ssflux',
        'latent','sensible','lwrad','swrad',
        'sustr','svstr','bustr','bvstr']

# check that all variables are actually there
for vn in vn_list:
    if vn not in vnl0:
        if Ldir['verbose']:
            print('ERROR: ds is missing %s, removing from vn_list')
        vn_list.remove(vn)

# override for testing
if Ldir['testing']:
    if Ldir['verbose']:
        print('Testing')
        print('- variables that would have been processed:')
        print(vn_list)
    vn_list = ['zeta','salt','temp','u','v']

for ii in range(Ldir['ii0'], Ldir['ii1']+1):
    ds = xr.open_dataset(fn_list[ii])
    if ii == Ldir['ii0']:
        lp = (ds[vn_list] * gs[ii]).squeeze(dim='ocean_time', drop=True).compute()
    else:
        lp = (lp + (ds[vn_list] * gs[ii]).squeeze(dim='ocean_time', drop=True)).compute()

# save the temp file to NetCDF
fstr = ('000' + str(Ldir['fnum']))[-2:]
lp.to_netcdf(temp_out_dir / ('lp_temp_' + fstr + '.nc'))

