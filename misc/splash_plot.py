"""
A pretty plot for talks.
"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
import Lfun
import zrfun
import zfun
import numpy as np

import plotting_functions as pfun
import matplotlib.pyplot as plt
import xarray as xr

Ldir = Lfun.Lstart()

# model output
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2019.07.04' / 'ocean_his_0001.nc')
xs = xr.open_dataset(fn)
x = xs.lon_psi.values
y = xs.lat_psi.values
th = xs['temp'][0,-1,1:-1,1:-1].values

# topography
tfn = (Ldir['parent'] / 'ptools_data' / 'topo' /
    'srtm15' / 'topo15.nc')
txs = xr.open_dataset(tfn)
step = 3
tx = txs['lon'][::step].values
ty = txs['lat'][::step].values
tz = txs['z'][::step,::step].values
tz[tz<0] = np.nan

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(8,12))
fig = plt.figure()

ax = fig.add_subplot(111)

cmap = 'RdYlBu_r'
cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
#fig.colorbar(cs)
    
cmap = 'gist_earth'
cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
#fig.colorbar(cs)


pfun.add_coast(ax)
pfun.dar(ax)
ax.axis([-130, -122, 42, 52])
ax.set_axis_off()
fig.tight_layout()

out_dir = Ldir['LOo'] / 'misc'
Lfun.make_dir(out_dir)
plt.savefig(out_dir / 'splash.png', transparent=True)

plt.show()
pfun.end_plot()