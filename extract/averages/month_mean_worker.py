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
parser.add_argument('-this_ym', type=str) # e.g. 2019_07
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-days_per_month', type=int) # number of days in this month
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
temp_out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'averages' / ('temp_' + Ldir['this_ym'])
# Lfun.make_dir(temp_out_dir) # Assume it was created by the calling function.

# set up to do the partial mean
fn_list = Lfun.get_fn_list('lowpass', Ldir, Ldir['ds0'], Ldir['ds1'])
if Ldir['verbose']:
    print(fn_list[0])
    print(fn_list[-1])

averaging_factor = 1/Ldir['days_per_month']

ii = 0
for fn in fn_list:
    ds = xr.open_dataset(fn)
    if ii == 0:
        lp = (ds*averaging_factor).compute()
    else:
        lp = (lp + ds*averaging_factor).compute()
    ii += 1

# save the temp file to NetCDF
fstr = ('000' + str(Ldir['fnum']))[-2:]
lp.to_netcdf(temp_out_dir / ('mean_temp_' + fstr + '.nc'))

