"""
Generic code to plot any mooring extraction
"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import plotting_functions as pfun

import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
moor_name = Lfun.choose_item(in_dir, tag='.nc', exclude_tag='', itext='** Choose mooring extraction from list **')
moor_fn = in_dir / moor_name

ds = nc.Dataset(moor_fn)

ot = ds['ocean_time'][:]
print('time step of mooring'.center(60,'-'))
print(np.diff(ot))

tind = [Lfun.modtime_to_datetime(tt) for tt in ot]
print('time limits'.center(60,'-'))
print('start ' + str(tind[0]))
print('end   ' + str(tind[-1]))

print('info'.center(60,'-'))
VN_list = []
for vn in ds.variables:
    print('%s %s' % (vn, ds[vn].shape))
    VN_list.append(vn)
    
# populate list of variables to plot
vn_list = []
if 'salt' in VN_list:
    vn_list += ['salt', 'temp']
if 'NO3' in VN_list:
    vn_list += ['NO3', 'oxygen', 'phytoplankton']

df = pd.DataFrame(index=tind)
df['zeta'] = ds['zeta'][:,0,0]
for vn in vn_list:
    df[vn] = ds[vn][:, -1, 0, 0]

pfun.start_plot()

df.plot(subplots=True, figsize=(16,10))

pfun.end_plot()