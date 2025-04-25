"""
Code to download a sequence of GLORYS hindcast daily averages.

For region7 each output file is 7 MB and takes about 1.5 minutes on my mac,
and three minutes on apogee.

We send the output to LO_data/glorys to match the pattern we used for hycom.
"""

import sys
import argparse
from datetime import datetime, timedelta
from lo_tools import Lfun
from lo_tools import glorys_functions as gfun
from time import time

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-src', '--source', type=str) # hindcast, interim, or forecast
parser.add_argument('-r', '--region', type=str, default='region7')
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
source = args.source
region = args.region

# check for required arguments
argsd = args.__dict__
for a in ['ds0','source']:
    if argsd[a] == None:
        print('*** Missing required argument for get_glorys_days.py: ' + a)
        sys.exit()

verbose = False
if args.testing:
    verbose = True
    from importlib import reload
    reload(gfun)

# get Ldir
Ldir = Lfun.Lstart()

# set time range to process
ds0 = args.ds0
if len(args.ds1) == 0:
    ds1 = ds0
else:
    ds1 = args.ds1
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

# loop over all days
dt = dt0
while dt <= dt1:
    tt0 = time()
    dstr = dt.strftime(Lfun.ds_fmt)
    if source in ['hindcast','interim']:
        # creates a single file
        out_dir = Ldir['data'] / 'glorys' / region
        Lfun.make_dir(out_dir)
        gfun.get_glorys_hindcast(dt, out_dir, source, region, verbose=verbose)
    elif source == 'forecast':
        out_dir = Ldir['data'] / 'glorys' / region / ('forecast_' + dstr)
        Lfun.make_dir(out_dir)
        gfun.get_glorys_forecast(dt, out_dir, region, verbose=verbose)
    dt += timedelta(days=1)
    print('Time for download = %0.1f\n' % (time()-tt0))
    sys.stdout.flush()