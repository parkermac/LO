"""
Plot a map of Froude Number.  Potentially this is a piece of evidence to use for
understanding Salish Sea exchange flow.

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
from cmocean import cm

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
Ldir['roms_out'] = Ldir['roms_out2']

in_dir = Ldir['roms_out'] / Ldir['gtagex'] / 'f2019.07.04'
fn = in_dir / 'ocean_his_0001.nc'

ds = xr.open_dataset(fn)

# PLOT CODE
xrho = ds['lon_rho'][0,:].values
yrho = ds['lat_rho'][:,0].values

# define box
aa = [-123.25, -122.1, 47, 48.75]
ix0 = zfun.find_nearest_ind(xrho, aa[0])
ix1 = zfun.find_nearest_ind(xrho, aa[1])
iy0 = zfun.find_nearest_ind(yrho, aa[2])
iy1 = zfun.find_nearest_ind(yrho, aa[3])

h = ds.h[iy0:iy1, ix0:ix1].values
rho_bot = ds.rho[0, 0, iy0:iy1, ix0:ix1].values
rho_top = ds.rho[0, -1, iy0:iy1, ix0:ix1].values
drho = rho_bot - rho_top
u = ds.ubar[0, iy0:iy1, ix0-1:ix1].values
v = ds.vbar[0, iy0-1:iy1, ix0:ix1].values
u[np.isnan(u)] = 0
v[np.isnan(v)] = 0
uu = (u[:, 1:] + u[:, :-1])/2
vv = (v[1:, :] + v[:-1, :])/2
spd2 = uu**2 + vv**2
spd2[np.isnan(drho)] = np.nan
spd2[spd2 < .001] = .001 # avoid divide by zero errors

# approximate Richardson number
rho0 = ds.rho0.values
g = 9.8
Ri = g * drho * h / (rho0 * spd2)

# psi_grid coordinates
x, y = np.meshgrid(ds.lon_u.values[0,ix0-1:ix1], ds.lat_v.values[iy0-1:iy1,0])

# PLOTTING
plt.close('all')
pfun.start_plot(fs=10, figsize=(20,10))
fig = plt.figure()

xt = [-123.2, -122.2]
yt = [47, 47.5, 48, 48.5]

ax = fig.add_subplot(131)
cs = ax.pcolormesh(x, y, drho, vmin=0, vmax=5, cmap=cm.dense)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title(r'$\Delta\rho\ [kg\ m^{-3}]$')
ax.set_xticks(xt)
ax.set_yticks(yt)

ax = fig.add_subplot(132)
cs = ax.pcolormesh(x, y, np.sqrt(spd2), vmin=0, vmax=2, cmap=cm.speed)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title(r'Speed $[m\ s^{-1}]$')
ax.set_xticks(xt)
ax.set_yticks(yt)
ax.set_yticklabels([])

ax = fig.add_subplot(133)
cs = ax.pcolormesh(x, y, 4*Ri, vmin=0, vmax = 2, cmap='RdYlBu_r')
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title(r'$4 x Ri$')
ax.set_xticks(xt)
ax.set_yticks(yt)
ax.set_yticklabels([])

plt.show()
pfun.end_plot()