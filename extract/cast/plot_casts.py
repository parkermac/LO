"""
Code to plot the results of a cast extraction.

run plot_casts -g cas6 -t v3 -x lo8b -ro 2 -cruises newport_line
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

import netCDF4 as nc
import Lfun
import plotting_functions as pfun
import matplotlib.pyplot as plt

in_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'cast' / Ldir['cruises']

fn_list = list(in_dir.glob('*.nc'))

foo = nc.Dataset(fn_list[0])
for vn in foo.variables:
    print('%14s: %s' % (vn, str(foo[vn].shape)))
foo.close()

plt.close('all')
pfun.start_plot(fs=14, figsize=(14,10))
fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False)
for fn in fn_list:
    ds = nc.Dataset(fn)
    axes[0,0].plot(ds['lon_rho'][0,0], ds['lat_rho'][0,0],'ob')
    axes[0,1].plot(ds['salt'][0,:,0,0], ds['temp'][0,:,0,0],'.g')
pfun.add_coast(axes[0,0])
pfun.dar(axes[0,0])
axes[0,1].set_xlabel('Salinity')
axes[0,1].set_ylabel('Potential Temperature')
axes[0,0].set_title(in_dir.name)
plt.show()
pfun.end_plot()
