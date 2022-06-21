"""
Generic code to plot any mooring extraction
"""
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
# you can choose either a file or a directory
moor_name = Lfun.choose_item(in_dir, itext='** Choose mooring extraction or folder from list **')
moor_item = in_dir / moor_name
if moor_item.is_file() and moor_name[-3:]=='.nc':
    moor_fn = moor_item
elif moor_item.is_dir():
    moor_name = Lfun.choose_item(moor_item, tag='.nc', itext='** Choose mooring extraction from list **')
    moor_fn = moor_item / moor_name

# load everything using xarray
ds = xr.open_dataset(moor_fn)
ot = ds.ocean_time.values
ot_dt = pd.to_datetime(ot)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
T = t/86400 # time in days from start
print('time step of mooring'.center(60,'-'))
print(t[1])
print('time limits'.center(60,'-'))
print('start ' + str(ot_dt[0]))
print('end   ' + str(ot_dt[-1]))
print('info'.center(60,'-'))
VN_list = []
for vn in ds.data_vars:
    print('%s %s' % (vn, ds[vn].shape))
    VN_list.append(vn)
    
# populate lists of variables to plot
vn2_list = ['zeta']
if 'Pair' in VN_list:
    vn2_list += ['shflux', 'swrad']
if 'ubar' in VN_list:
    vn2_list += ['ubar', 'vbar']
vn3_list = []
if 'salt' in VN_list:
    vn3_list += ['salt', 'temp']
if 'NO3' in VN_list:
    vn3_list += ['oxygen','NO3']
# if 'u' in VN_list:
#     vn3_list += ['u', 'v']

# drop missing variables
vn2_list = [item for item in vn2_list if item in ds.data_vars]
vn3_list = [item for item in vn3_list if item in ds.data_vars]

# plot time series using a pandas DataFrame
df = pd.DataFrame(index=ot)
for vn in vn2_list:
    df[vn] = ds[vn].values
for vn in vn3_list:
    df[vn] = ds[vn][:, -1]
    #df[vn] = zfun.lowpass(ds[vn][:, -1].values, f='godin')

plt.close('all')
pfun.start_plot()
df.plot(subplots=True, figsize=(16,12), grid=True)
plt.show()
pfun.end_plot()