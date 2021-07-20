"""
This is a one-off program to convert the output of a segment extraction
from an xarray DataArray to a Dataset.
"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))

import Lfun
import xarray as xr

Ldir = Lfun.Lstart(gridname='cas6', tag='v3', ex_name='lo8b')
year_str = str(2018)
date_str = '_' + year_str + '.01.01_' + year_str + '.12.31'
seg_fn = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / ('segments' + date_str + '.nc')
seg_ds_fn = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / ('segments_ds' + date_str + '.nc')
da = xr.open_dataarray(seg_fn)

vns = da.coords['vn'].values
segs = da.coords['seg'].values
times = da.coords['time'].values

ds = xr.Dataset(coords={'time': times,'seg': segs})

for vn in vns:
    v = da.sel(vn=vn).values
    ds[vn] = (['time','seg'], v)
    
ds.to_netcdf(seg_ds_fn)
