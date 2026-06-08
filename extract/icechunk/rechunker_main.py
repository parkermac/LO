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

run rechunker_main -gtx cas7_t2_x11b -0 2024.01.01 -1 2024.01.02 -lt hourly -test True

run rechunker_main -gtx cas7_t1_x11ab -0 2024.01.01 -1 2024.02.29 -lt average -test True

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

# pid = os.getpid()
# print(' rechunker_main '.center(60,'='))
# print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t2_x11b
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str)# hourly, average
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

# create input file list and subsets

if args.list_type == 'hourly':
    time_bin_list = pd.date_range(args.ds0,args.ds1,freq='D')
    for dt in time_bin_list:
        dt_list = pd.date_range(dt,dt+timedelta(days=1),freq='h')
        print('')
        print(dt_list)
elif args.list_type == 'average':
    time_bin_list = pd.date_range(args.ds0,args.ds1,freq='MS')
    me_list = pd.date_range(args.ds0,args.ds1,freq='ME')
    bin_dict = dict(zip(time_bin_list, me_list))
    for dt in time_bin_list:
        dt_list = pd.date_range(dt,bin_dict[dt],freq='D')
        print('')
        print(dt_list)


#print(time_bin_list)

#fn_list0 = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'], his_num=1)

# if Ldir['list_type'] == 'hourly':

# elif Ldir['list_type'] == 'average':

# else:
#     print('unsupported list type')
#     sys.exit()

