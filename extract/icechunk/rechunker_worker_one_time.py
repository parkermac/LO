"""
Process a single history file.
"""

import os
import xarray as xr
import pandas as pd
from time import time
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import argparse
import sys
from lo_tools import Lfun

# command line arguments
"""
-tid ${SLURM_ARRAY_TASK_ID}
indir
outdir
gtagex
this_ds0
this_ds1
list_type
"""
parser = argparse.ArgumentParser()

parser.add_argument('-tid', type=str)
parser.add_argument('-indir', type=str)
parser.add_argument('-outdir', type=str)
parser.add_argument('-gtagex', type=str)
parser.add_argument('-this_ds0', type=str)
parser.add_argument('-this_ds1', type=str)
parser.add_argument('-list_type', type=str)

args = parser.parse_args()

Ldir = Lfun.Lstart()

# 1. Figure out which file to get and copy it to /var/tmp using s5cmd

list_type = args.list_type
dt0 = datetime.strptime(args.this_ds0, Ldir['ds_fmt'])
dt1 = datetime.strptime(args.this_ds1, Ldir['ds_fmt'])

# zero-based task ID to use for indexing
tid0 = int(args.tid) - 1

# Associate a history file with a task ID
if list_type == 'hourly0':
    sub_ind = pd.date_range(dt0,dt1,freq='h')
elif list_type == 'hourly':
    sub_ind = pd.date_range(dt0+timedelta(hours=1),dt1,freq='h')
elif list_type == 'average':
    sub_ind = pd.date_range(dt0,dt1,freq='D')
ntimes = len(sub_ind)

this_dt = sub_ind[tid0]

if list_type in ['hourly0', 'hourly']:
    f_string = 'f' + this_dt.strftime(Ldir['ds_fmt'])
    hour = this_dt.hour
    nhis = ('000' + str(int(hour)+1))[-4:]
    fn0 = 'ocean_his_' + nhis + '.nc'
elif list_type == 'average':
    f_string = 'f' + this_dt.strftime(Ldir['ds_fmt'])
    fn0 = 'ocean_avg_0001.nc'
fn = args.indir + '/' + f_string + '/' + fn0

cmd_list = ['s5cmd','cp', fn, os.environ['TMPDIR'] + '/']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()

# 2. Use xarray or ncks to make a new Dataset with only selected variables 
# and save it to a folder in /gscratch

get_tsa = True
get_vel = True
get_bio = True
get_surfbot = True

# create a variable list focused on things we could usually compare to observations

# bio_list = ',NO3,NH4,phytoplankton,zooplankton,SdetritusN,LdetritusN,SdetritusC,LdetritusC,oxygen,alkalinity,TIC,rho'
bio_list = ',NO3,NH4,phytoplankton,oxygen,alkalinity,TIC'
# check to see if the model has WETDRY in use
do_wetdry = False
# if 'wetdry_mask_rho' in ds.data_vars:
#     do_wetdry = True

vn_list = 'h,zeta'
if do_wetdry:
    vn_list += ',wetdry_mask_rho'
if get_tsa:
    # vn_list += ',salt,temp,AKs,AKv'
    vn_list += ',salt,temp'
if get_vel:
    # vn_list += ',u,v,w,ubar,vbar'
    vn_list += ',u,v'
    if do_wetdry:
        vn_list += ',wetdry_mask_u,wetdry_mask_v'
if get_bio:
    vn_list += bio_list
if get_surfbot:
    vn_list += ',Pair,Uwind,Vwind,shflux,ssflux,latent,sensible,lwrad,swrad,sustr,svstr,bustr,bvstr'

in_fn = os.environ['TMPDIR'] + '/' + fn0
padded_tid = ('000' + args.tid)[-4:]
out_fn = args.outdir + '/subset_' + padded_tid + '.nc'

if False:
    # method 1: use xarray
    ds = xr.open_dataset(in_fn)
    # do a final check to drop missing variables from the list
    vn_list = (',').join([item for item in vn_list.split(',') if item in ds.data_vars])
    ds_out = ds[vn_list.split(',')]
    ds_out.to_netcdf(out_fn)
else:
    # Method 2: use ncks
    cmd_list = ['ncks',
            '-v', vn_list,
            '--mk_rec_dim', 'ocean_time','-O', str(in_fn), str(out_fn)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc.communicate()


