"""

Extract as-run river time series for the cas6_v0_live forecast, which exists
for the time range 2016.12.15 to 2023.05.21.

Includes all NPZD tracers, and package the results as
an xarray Dataset.

Output: LO_output/pre/river1/lo_base/extraction_2016.12.15_2023.05.21.nc

This code ONLY works for this run, and assumes you are running on apogee, so
we skip command line arguments and hard-code choices about the run.

To run on apogee:
run extract_rivers_cas6_v0_live

Performance: takes 25 sec per year on apogee.

NOTE: The forcing files were originally spread across two machines and different naming
conventions with the change happening around October 2021.
Because one of the machines, boiler, is being decommissioned, on 2023.09.25 I created
LO/misc/copy_river_files.py and ran it on apogee. This copied
all the rivers.nc files for 2016.12.15 to 2021.10.16
 - from /boildat/parker/LiveOcean_output/cas6_v3/f[]/riv2
 - to /dat1/parker/LO_output/cas6_v0/f[]/riv0
so now there is a continuous river forcing record on apogee for the cas6_v0_live run
for its entire length 2016.12.15 to 2023.05.21

"""

from lo_tools import Lfun, zrfun
from datetime import datetime, timedelta
from time import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
    
Ldir = Lfun.Lstart(gridname='cas6',tag='v0',ex_name='live')
Ldir['ds0'] = '2016.12.15'
Ldir['ds1'] = '2023.05.21'
ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt0 = time()

# long list of variables to extract
vn_list = ['transport', 'salt', 'temp', 'oxygen',
    'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
    'TIC', 'alkalinity']

# input directory
in_dir = Ldir['LOo'] / 'forcing' / 'cas6_v0'

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'pre' / 'river1' / 'lo_base' / 'Data_roms'
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

# Get list of river names.
# This is a bit titchy because of NetCDF 3 limitations on strings, forcing them
# to be arrays of characters. This was the case for the start of this time series.
mds = mds_list[0]
fn = in_dir / ('f' + mds) / 'riv0' / 'rivers.nc'
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
    
    this_dt = datetime.strptime(mds, Lfun.ds_fmt)
    if this_dt.day == 1 and this_dt.month == 1:
        print(' Year = %d' % (this_dt.year))
        
    fn = in_dir / ('f' + mds) / 'riv0' / 'rivers.nc'
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
    
