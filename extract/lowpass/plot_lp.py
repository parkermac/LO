"""
Code to make a plot that informally tests the lowpass code.

Result: looks good by eye.
"""

from lo_tools import Lfun, zrfun, zfun
import xarray as xr
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart(gridname='cas6', tag='v0', ex_name = 'live')

his_fn = Ldir['roms_out'] / Ldir['gtagex'] / 'f2019.07.05' / 'ocean_his_0013.nc'
lp_fn = Ldir['roms_out'] / Ldir['gtagex'] / 'f2019.07.05' / 'lowpassed.nc'

his_ds = xr.open_dataset(his_fn)
lp_ds = xr.open_dataset(lp_fn)

plon, plat = pfun.get_plon_plat(his_ds.lon_rho.values, his_ds.lat_rho.values)

vmin = 28
vmax = 33
cmap='Spectral_r'

plt.close('all')
pfun.start_plot(figsize=(18,12))
fig = plt.figure()

ax = fig.add_subplot(121)
cs = ax.pcolormesh(plon,plat,lp_ds.salt[0,-1,:,:].values, cmap=cmap, vmin=vmin, vmax=vmax)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(lp_ds))
ax.set_title('lowpassed')

ax = fig.add_subplot(122)
cs = ax.pcolormesh(plon,plat,his_ds.salt[0,-1,:,:].values, cmap=cmap, vmin=vmin, vmax=vmax)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(his_ds))
ax.set_title('history')

plt.show()
pfun.end_plot()
