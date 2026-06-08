"""
Code to go through a sequence of model history files, divide it into smaller subsets of files, and:
(i) select a subset of varibles
(ii) concatenate a bunch of times into a single NetCDF file
(iii) add z variables
(iv) rechunk the result so it has a larger chunksize for ocean_time
(v) save to kopah for later processing using icechunk

This is designed to work on klone, accessing model output from kopah.

The smaller subset will be all the hours in a day for hourly list_type, or all the days in a month for average list_type.

Testing:

run rechunker_main -gtx cas7_t2_x11b -0 2024.01.01 -1 2024.01.02 -lt hourly0 -test True

run rechunker_main -gtx cas7_t1_x11ab -0 2024.01.01 -1 2024.02.29 -lt average -test True

Note, using list_type = hourly0 zero means that on the first day we will get 25 hours (0:24)
and on subsequent days it will revert to hourly, getting 24 hours (1:24).

When populating a NEW collection you should use hourly0
but when adding to it later you should use hourly.

"""

# imports
import sys
import argparse
from lo_tools import Lfun, zfun, zrfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
from time import time
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta

pid = os.getpid()
print(' rechunker_main '.center(60,'='))
print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t2_x11b
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str)# hourly, hourly0, average
# set this to True to interpolate all u, and v fields to the rho-grid
parser.add_argument('-uv_to_rho', default=True, type=Lfun.boolean_string)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()

Ldir = Lfun.Lstart()

if args.testing:
    Ldir['local_user'] = 'pmacc'

# name the s3 bucket to look in
bucket_name = 'liveocean-' + Ldir['local_user']

# create start and end day indexes where each item is a start or end day
# for a given subset

if args.list_type in ['hourly','hourly0']:
    dt0_ind = pd.date_range(args.ds0,args.ds1,freq='D') # a DatetimeIndex
    dt1_ind = dt0_ind + timedelta(days=1)
elif args.list_type == 'average':
    dt0_ind = pd.date_range(args.ds0,args.ds1,freq='MS')
    dt1_ind = pd.date_range(args.ds0,args.ds1,freq='ME')
dt_dict = dict(zip(dt0_ind, dt1_ind))

# name some things to pass to the subprocess
gtagex = args.gtagex
indir = 's3://' + bucket_name + '/LO_roms/' + gtagex
outdir = str(Ldir['LOo']) + '/icechunk_temp'
list_type = args.list_type

# create a clean directory for the output
if not args.testing:
    Lfun.make_dir(outdir, clean=True)

# loop over the subsets
for dt0 in dt0_ind:

    if (list_type == 'hourly0') and (dt0 > dt0_ind[0]):
        list_type = 'hourly'

    dt1 = dt_dict[dt0]

    # find the number of jobs to run per subset
    if list_type == 'hourly0':
        sub_ind = pd.date_range(dt0,dt1,freq='h')
    elif list_type == 'hourly':
        sub_ind = pd.date_range(dt0+timedelta(hours=1),dt1,freq='h')
    elif list_type == 'average':
        sub_ind = pd.date_range(dt0,dt1,freq='D')
    ntimes = len(sub_ind)

    this_ds0 = dt0.strftime(Ldir['ds_fmt'])
    this_ds1 = dt1.strftime(Ldir['ds_fmt'])
    cmd_list = ['sbatch','--array=1-'+str(ntimes),
                './rechunker_worker.sh',
                indir, outdir, gtagex,
                this_ds0, this_ds1, list_type]
    
    if args.testing:
        print('')
        print(' '.join(cmd_list))
        print(sub_ind)
    else:
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\nSTDOUT')
            print(stdout.decode())
        if len(stderr) > 0:
            print('\nSTDERR')
            print(stderr.decode())

# dt_list = pd.date_range(dt,dt+timedelta(days=1),freq='h')
# dt_list = pd.date_range(dt,bin_dict[dt],freq='D')

