"""

Extract as-run river time series.

To run on boiler:
run extract_rivers -g cas6 -t v3 -x junk -0 2018.01.01 -1 2018.01.10
run extract_rivers -g cas6 -t v3 -x junk -0 2018.01.01 -1 2018.12.31

Performance: takes 3 seconds per year on boiler

Modified to include all NPZD tracers, and package the results as
an xarray Dataset.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from time import time
import Lfun
import numpy as np
import netCDF4 as nc
import pandas as pd
import zrfun
import xarray as xr
    
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
fn = Ldir['parent'] / 'LiveOcean_output' / Ldir['gtag'] / ('f' + mds) / 'riv2' / 'rivers.nc'
ds = nc.Dataset(fn)
rn = ds['river_name'][:]
NR = rn.shape[1]
riv_name_list = []
for ii in range(NR):
    a = rn[:,ii]
    if isinstance(a, np.ma.MaskedArray):
        a = a.data
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
    v_dict[vn] = nanmat
tt = 0
for mds in mds_list:
    fn = Ldir['parent'] / 'LiveOcean_output' / Ldir['gtag'] / ('f' + mds) / 'riv2' / 'rivers.nc'
    ds = nc.Dataset(fn)
    # The river transport is given at noon of a number of days surrounding the forcing date.
    # Here we find the index of the time for the day "mds".
    RT = ds['river_time'][:]
    ii = 0
    for rt in RT:
        rdt = Lfun.modtime_to_datetime(rt)
        rds = datetime.strftime(rdt, Lfun.ds_fmt)
        if rds == mds:
            for vn in vn_list:
            v_dict[vn][tt,:] = ds['river_' + vn][ii,:]
        ii += 1
    ds.close()
    tt += 1

# make transport positive    
v_dict['river_transport'] = np.abs(v_dict['river_transport'])

# store output in an xarray Dataset
mdt_list = [(datetime.strptime(item, Lfun.ds_fmt) + timedelta(days=0.5)) for item in mds_list]
times = pd.Index(mdt_list)

ds = xr.Dataset(coords={'time': times,'riv': river_name_list})

for vn in vn_list:
    v = v_dict[vn]
    ds[vn] = (['time','riv'], v)
    
ds.to_netcdf(out_fn)

print('Total time for extraction = %d seconds' % (time() - tt0))
    
