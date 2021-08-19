"""
Code to plot the results of a cast extraction.

run plot_casts -gtx cas6_v3_lo8b -cruises newport_line
"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

import xarray as xr
import matplotlib.pyplot as plt
from time import time

in_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'cast' / Ldir['cruises']

fn_list = list(in_dir.glob('*.nc'))

foo = xr.open_dataset(fn_list[0])
for vn in foo.data_vars:
    print('%14s: %s' % (vn, str(foo[vn].shape)))
foo.close()

tt0 = time()
x = []; y = []
s0 = []; s1 = []
t0 = []; t1 = []
for fn in fn_list:
    ds = xr.open_dataset(fn)
    x.append(ds.lon_rho.values)
    y.append(ds.lat_rho.values)
    s0.append(ds.salt[0].values)
    t0.append(ds.temp[0].values)
    s1.append(ds.salt[-1].values)
    t1.append(ds.temp[-1].values)
    ds.close()
print('Took %0.2f sec' % (time()-tt0))

plt.close('all')
pfun.start_plot(fs=14, figsize=(14,10))
fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False)
axes[0,0].plot(x,y,'og')
axes[0,1].plot(s0,t0,'.', c='orange')
axes[0,1].plot(s1,t1,'.', c='dodgerblue')
pfun.add_coast(axes[0,0])
pfun.dar(axes[0,0])
axes[0,1].set_xlabel('Salinity')
axes[0,1].set_ylabel('Potential Temperature')
axes[0,0].set_title(in_dir.name)
plt.show()
pfun.end_plot()
