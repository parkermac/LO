"""
This code is to report on the model forcing. It prints out the
max and min values of all variables in all forcing files that
have time in one of their dimensions.

I used this to diagnose a forcing problem that was keeping the
long hindcast from running past 2015.09.23. It revealed that the
min value for Tair in the atm00 forcing for cas7 was -273.15 degrees
which is pretty cold!

Note

Example call:
run forcing_report -g cas7 -0 2015.09.23

"""

from lo_tools import Lfun
import xarray as xr
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2015.09.23
args = parser.parse_args()
gridname = args.gridname
fstr = 'f' + args.ds0

Ldir = Lfun.Lstart()
in_dir = Ldir['LOo'] / 'forcing' / gridname / fstr

d_list = [f for f in in_dir.iterdir() if f.is_dir()]

for d in d_list:
    print('\n%s' % (str(d)))
    nc_list = d.glob('*.nc')
    for fn in nc_list:
        print('  %s' % (fn.name))
        ds = xr.open_dataset(fn)
        vlist = [v for v in ds.data_vars]
        for v in vlist:
            dim_list = [dim for dim in ds[v].dims if 'time' in dim]
            if len(dim_list) > 0:
                print('    %s %s' % (v, ds[v].dims))
                print('     max=%0.3f min=%0.3f' % (ds[v].max(), ds[v].min()))
        ds.close()
