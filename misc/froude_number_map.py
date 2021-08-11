"""
Plot a map of Froude Number.  Potentially this is a piece of evidence to use for
understanding Salish Sea exchange flow.

"""

from pathlib import Path
import sys

pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr

import plotting_functions as pfun
import zfun
import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
Ldir['roms_out'] = Ldir['roms_out2']

in_dir = Ldir['roms_out'] / Ldir['gtagex'] / 'f2019.07.04'
fn = in_dir / 'ocean_his_0001.nc'

ds = xr.open_dataset(fn)
x = ds.lon_psi.values
y = ds.lat_psi.values

ix0 = zfun.find_nearest_ind(x[0,:], -125.5)
ix1 = zfun.find_nearest_ind(x[0,:], -122)
iy0 = zfun.find_nearest_ind(y[:,0], 47)
iy1 = zfun.find_nearest_ind(y[:,0], 50)

x = x[iy0:iy1+1,ix0:ix1+1]
y = y[iy0:iy1+1,ix0:ix1+1]

DR = (ds.rho[0,0,:,:] - ds.rho[0,-1,:,:]).values
DR = DR[iy0+1:iy1+1,ix0+1:ix1+1]

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(20,10))
fig = plt.figure()

ax = fig.add_subplot(121)
cs = ax.pcolormesh(x, y, DR)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)

plt.show()
pfun.end_plot()