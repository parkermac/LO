"""
Code to find a time series of the time step for a sequence of days.

Testing on mac:

run timestep_record.py

To run for real on perigee
python timestep_record.py -gtx cas6_v3_lo8b -ro 2 -0 2017.01.01 -1 2021.11.30

This runs very quickly, like a few seconds.
"""

import argparse
from datetime import datetime, timedelta
import pandas as pd

from lo_tools import Lfun

parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str, default='cas6_v3_lo8b')
parser.add_argument('-ro', '--roms_out_num', type=int, default=2)
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str, default='2019.07.04')
parser.add_argument('-1', '--ds1', type=str, default='2019.07.06')

# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__

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
    
ds0 = args.ds0
if len(args.ds1) == 0:
    ds1 = ds0
else:
    ds1 = args.ds1
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

# Loop over days
dt_ser = pd.Series(dtype=float)
dt = dt0
ii = 0
while dt <= dt1:
    f_string = 'f' + dt.strftime(Lfun.ds_fmt)
    roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    fn = roms_out_dir / 'log.txt'
    # example line I am searching for
    #      40.000  dt                Timestep size (s) for 3-D equations.
    try:
        with open(fn, 'r') as x:
            for line in x:
                if 'Timestep size (s) for 3-D equations' in line:
                    timestep = line.strip().split(' ')[0]
                    dt_ser[dt] = float(timestep)
                    ii += 1
                    break
    except FileNotFoundError:
        pass
    dt += timedelta(days=1)

print(ii)

# save_output
out_dir = Ldir['LOo'] / 'misc'
Lfun.make_dir(out_dir)
out_fn = out_dir / ('dt_ser_' + Ldir['gtagex'] + '.p')
out_fn.unlink(missing_ok=True)
dt_ser.to_pickle(out_fn)

print('Pickled Series saved in:')
print(str(out_fn))

