"""

Extract as-run river time series.

To test on mac:
run extract_rivers -gtx cas6_v3_lo8b -0 2019.07.04 -1 2019.07.04

To run on boiler:
run extract_rivers -gtx cas6_v3_lo8b -0 2018.01.01 -1 2018.01.10
run extract_rivers -gtx cas6_v3_lo8b -0 2018.01.01 -1 2018.12.31

Performance: takes 23 sec per year on boiler

Modified to include all NPZD tracers, and package the results as
an xarray Dataset.

***
NOTE: this is hard-coded to LiveOcean_output / [gtag] / riv2 so it
pretty specific to the cas6_v3_lo8b run.  Also, it expects to find all
the NPZDOC variables.
***

"""

from lo_tools import Lfun, zrfun
from lo_tools import extract_argfun as exfun
    
Ldir = exfun.intro() # this handles the argument passing

from datetime import datetime, timedelta
from time import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
    
ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt0 = time()

# long list of variables to extract
vn_list = ['transport', 'salt', 'temp', 'oxygen',
    'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
    'TIC', 'alkalinity']

print(' Doing river extraction for '.center(60,'='))
print(' gtag = ' + Ldir['gtag'])
outname = 'extraction_' + ds0 + '_' + ds1 + '.nc'

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms'
Lfun.make_dir(out_dir)

out_fn = out_dir / outname
out_fn.unlink(missing_ok=True)

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
ndays = (dt1-dt0).days + 1

# make mds_list: list of datestrings (e.g. 2017.01.01) to loop over
mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, Lfun.ds_fmt))
    mdt = mdt + timedelta(days=1)

# get list of river names
# (this is a bit titchy because of NetCDF 3 limitations on strings, forcing them
# to be arrays of characters)
mds = mds_list[0]
fn = Path('/boildat1').absolute() / 'parker' / 'LiveOcean_output' / Ldir['gtag'] / ('f' + mds) / 'riv2' / 'rivers.nc'
ds = xr.open_dataset(fn)
rn = ds['river_name'].values
NR = rn.shape[1]
riv_name_list = []
for ii in range(NR):
    a = rn[:,ii]
    r = []
    for l in a:
        r.append(l.decode())
    rr = ''.join(r)
    riv_name_list.append(rr)
ds.close()

NT = len(mds_list)

nanmat = np.nan * np.ones((NT, NR))
v_dict = dict()
for vn in vn_list:
    v_dict[vn] = nanmat.copy()
tt = 0
for mds in mds_list:
    fn = Path('/boildat1').absolute() / 'parker' / 'LiveOcean_output' / Ldir['gtag'] / ('f' + mds) / 'riv2' / 'rivers.nc'
    ds = xr.open_dataset(fn)
    # The river transport is given at noon of a number of days surrounding the forcing date.
    # Here we find the index of the time for the day "mds".
    RT = pd.to_datetime(ds['river_time'].values)
    mdt = datetime.strptime(mds, Lfun.ds_fmt) + timedelta(hours=12)
    mask = RT == mdt
    for vn in vn_list:
        if vn == 'transport':
            v_dict[vn][tt,:] = ds['river_' + vn][mask,:]
        else:
            # the rest of the variables allow for depth variation, but we
            # don't use this, so, just use the bottom value
            v_dict[vn][tt,:] = ds['river_' + vn][mask,0,:]
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
    
