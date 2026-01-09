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
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2015.09.23
parser.add_argument('-f', '--frc', type=str, default='all')
# You can use the -f flag to select a single forcing type. The default
# 'all' will loop over every forcing type in the day folder.
parser.add_argument('-l', '--list_only', type=Lfun.boolean_string, default=False)
# Use -l True to just list the forcing types available for this day.
parser.add_argument('-nan', '--show_nansum', type=Lfun.boolean_string, default=False)
# Use -nan True to include nansum values for each time step.
args = parser.parse_args()
gridname = args.gridname
fstr = 'f' + args.ds0
frc = args.frc

Ldir = Lfun.Lstart()
in_dir = Ldir['LOo'] / 'forcing' / gridname / fstr

d0_list = [f for f in in_dir.iterdir() if f.is_dir()]

if args.list_only:
    f_list = [f.name for f in d0_list]
    print(f_list)
    sys.exit()

if frc != 'all':
    d_list = [f for f in d0_list if f.name==frc]
elif frc == 'all':
    d_list = d0_list

if len(d_list) == 0:
    print('Requested forcing type: %s not present this day.' % (frc))
    print(d_list)
    sys.exit()

for d in d_list:
    print('\n%s' % (str(d)))
    nc_list = d.glob('*.nc')
    for fn in nc_list:
        if fn.name == 'ocean_bry.nc':
            # Skip this because it clutters up the output, and any errors
            # should be evident in the ocean_clm output.
            pass
        else:
            print('\n  %s' % (fn.name))
            ds = xr.open_dataset(fn, decode_times=False)
            vlist = [v for v in ds.data_vars]
            for v in vlist:
                dim_list = [dim for dim in ds[v].dims if 'time' in dim]
                if len(dim_list) > 0:
                    print('\n    %s %s' % (v, ds[v].dims))
                    print('     max=%0.3f min=%0.3f' % (ds[v].max(), ds[v].min()))
                    if args.show_nansum:
                        # assumes time is the first index
                        vv = ds[v].values
                        nt = vv.shape[0]
                        for tt in range(nt):
                            print('     - time %3d: %d nans' % (tt, int(np.sum(np.isnan(vv)))))

            ds.close()
