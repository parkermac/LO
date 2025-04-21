"""
Code to download a sequence of GLORYS hindcast daily averages.
"""

import sys
import argparse
from datetime import datetime, timedelta
from lo_tools import Lfun
from lo_tools import glorys_functions as gfun

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-r', '--region', type=str, default='region7')
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['ds0']:
    if argsd[a] == None:
        print('*** Missing required argument for get_glorys_days.py: ' + a)
        sys.exit()

if args.testing:
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

    dstr = dt.strftime(Lfun.ds_fmt)
    out_dir = Ldir['data'] / 'glorys' / args.region
    Lfun.make_dir(out_dir)
    gfun.get_glorys_hindcast(dt, out_dir, args.region, verbose=False)

    dt += timedelta(days=1)