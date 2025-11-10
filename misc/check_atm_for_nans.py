"""
Code to check atm forcing for nans. This is responding to the blowup problem
due to hour 19 of 2014.06.18 in atm00.

We will do the inspection by just looking for the word WARNING in the screen_output.
"""

import argparse
from datetime import datetime, timedelta
from lo_tools import Lfun
Ldir = Lfun.Lstart()

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-f', '--frc', type=str)        # e.g. atm00
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2014.08.17
parser.add_argument('-1', '--ds1', type=str)        # e.g. 2014.06.18
args = parser.parse_args()

dt0 = datetime.strptime(args.ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(args.ds1, Lfun.ds_fmt)

dt = dt0
while dt <= dt1:
    in_fn = Ldir['LOo'] / 'forcing' / args.gridname / ('f' + dt.strftime(Lfun.ds_fmt)) / args.frc / 'Info' / 'screen_output.txt'
    with open(in_fn, 'r') as file:
        for line in file:
            if 'WARNING' in line:
                print('Warning found in ' + dt.strftime(Lfun.ds_fmt))
                break

    dt += timedelta(days=1)
