"""
Utility code to check for missing days in a long hindcast.
"""

import pandas as pd
import argparse
from datetime import datetime
from lo_tools import Lfun

# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gtx', '--gtagex', type=str, default='cas7_t1_x11ab')
parser.add_argument('-ro', '--roms_out_num', type=int, default=1) # 1 = Ldir['roms_out1'], etc.
parser.add_argument('-hn', '--his_num', type=int, default=2)
parser.add_argument('-0', '--ds0', type=str, default='2012.10.07')
parser.add_argument('-1', '--ds1', type=str, default='2017.12.31')
args = parser.parse_args()

Ldir = Lfun.Lstart()
in_dir = Ldir['roms_out'] / args.gtagex

dt0 = datetime.strptime(args.ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(args.ds1, Lfun.ds_fmt)
his_str = ('0000' + str(args.his_num))[-4:]

dtr = pd.date_range(dt0,dt1)

print('\n' + str(in_dir))
print('Checking %s to %s\n' % (args.ds0, args.ds1))
for dt in dtr:
    f_str = 'f' + dt.strftime(Lfun.ds_fmt)
    fn = in_dir / f_str / ('ocean_his_' + his_str + '.nc')

    if not fn.is_file():
        print('Missing ' + f_str)