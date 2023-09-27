"""

Extract as-run river time series specific to a grid and a river forcing type.

Includes all NPZD tracers, and package the results as
an xarray Dataset.

It uses the gridname and the river forcing type to create a
"collection tag" = [ctag], e.g. cas6_riv00

Output: LO_output/pre/river1/[ctag]/extraction_[date range].nc

To test on mac:
run extract_rivers -g cas6 -0 2019.07.04 -1 2019.07.04 -riv riv00

To run on apogee:
run extract_rivers -g cas6 -0 2022.01.01 -1 2022.12.31 -riv riv00

Performance: 55 sec per year on apogee.

"""

from lo_tools import Lfun, zrfun
from datetime import datetime, timedelta
from time import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas6
parser.add_argument('-0', '--ds0', type=str) # e.g. 2022.01.01
parser.add_argument('-1', '--ds1', type=str) # e.g. 2022.12.31
parser.add_argument('-riv', type=str) # e.g. riv00
args = parser.parse_args()
argsd = args.__dict__
for a in ['gridname','ds0','ds1','riv']:
    if argsd[a] == None:
        print('*** Missing required argument to extract_argfun.intro(): ' + a)
        sys.exit()
Ldir = Lfun.Lstart(gridname=args.gridname)
ds0 = args.ds0
ds1 = args.ds1
ctag = Ldir['gridname'] + '_' + args.riv

tt0 = time()


# input directory
in_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname']

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag / 'Data_roms'
Lfun.make_dir(out_dir)
outname = 'extraction_' + ds0 + '_' + ds1 + '.nc'
out_fn = out_dir / outname
out_fn.unlink(missing_ok=True)

print(' Creating river extraction to file '.center(60,'='))
print(str(out_fn))

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
ndays = (dt1-dt0).days + 1

# make mds_list: list of datestrings (e.g. 2017.01.01) to loop over
mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, Lfun.ds_fmt))
    mdt = mdt + timedelta(days=1)

# Get list of river names and set variables to use.
#
# NOTE: This may not work with the old ROMS output, like the early years of
# cas6_v0_live. For that case use the custom code:
# extract_rivers_cas6_v0_live.py.
mds = mds_list[0]
fn = in_dir / ('f' + mds) / args.riv / 'rivers.nc'
ds = xr.open_dataset(fn)
riv_name_list = list(ds['river_name'].values)
NR = len(riv_name_list)
# long list of variables to extract
if 'river_NH4' in ds.data_vars:
    # new ROMS version
    vn_list = ['transport', 'salt', 'temp', 'Oxyg',
        'NH4','NO3', 'Phyt', 'Zoop', 'SDeN', 'LDeN',
        'TIC', 'TAlk']
else:
    # old ROMS version (like cas6_v0_live)
    # likely not needed, but retained to support legacy cases
    vn_list = ['transport', 'salt', 'temp', 'oxygen',
        'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
        'TIC', 'alkalinity']
ds.close()

NT = len(mds_list)

nanmat = np.nan * np.ones((NT, NR))
v_dict = dict()
for vn in vn_list:
    rvn = 'river_' + vn
    if rvn in ds.data_vars:
        v_dict[vn] = nanmat.copy()
    else:
        pass
tt = 0
for mds in mds_list:

    this_dt = datetime.strptime(mds, Lfun.ds_fmt)
    if this_dt.day == 1 and this_dt.month == 1:
        print(' Year = %d' % (this_dt.year))

    fn = in_dir / ('f' + mds) / args.riv / 'rivers.nc'
    ds = xr.open_dataset(fn)
    # The river transport is given at noon of a number of days surrounding the forcing date.
    # Here we find the index of the time for the day "mds".
    RT = pd.to_datetime(ds['river_time'].values)
    mdt = this_dt + timedelta(hours=12)
    mask = RT == mdt
    for vn in vn_list:
        if vn == 'transport':
            v_dict[vn][tt,:] = ds['river_' + vn][mask,:]
        else:
            rvn = 'river_' + vn
            if rvn in ds.data_vars:
                # the rest of the variables allow for depth variation, but we
                # don't use this, so, just use the bottom value
                v_dict[vn][tt,:] = ds['river_' + vn][mask,0,:]
            else:
                pass
    ds.close()
    tt += 1

# make transport positive
v_dict['transport'] = np.abs(v_dict['transport'])

# store output in an xarray Dataset
mdt_list = [(datetime.strptime(item, Lfun.ds_fmt) + timedelta(hours=12)) for item in mds_list]
times = pd.Index(mdt_list)

x = xr.Dataset(coords={'time': times,'riv': riv_name_list})

for vn in vn_list:
    v = v_dict[vn]
    x[vn] = (('time','riv'), v)

x.to_netcdf(out_fn)
x.close()

print('Total time for extraction = %d seconds' % (time() - tt0))
    
