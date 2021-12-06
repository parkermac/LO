"""
Code to compare the results of LO/forcing/tide0 to LiveOcean/forcing/tide2

Result - based on two spot checks they look identical.

"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn0 = '/Users/pm8/Documents/LO_output/forcing/cas6_v3/f2019.07.04/tide0/tides.nc'

fn2 = '/Users/pm8/Documents/LiveOcean_output/cas6_v3/f2019.07.04/tide2/tides.nc'

fng = '/Users/pm8/Documents/LO_data/grids/cas6/grid.nc'

ds0 = nc.Dataset(fn0)
ds2 = nc.Dataset(fn2)
dsg = nc.Dataset(fng)

x = dsg['lon_psi'][:]
y = dsg['lon_psi'][:]
m = dsg['mask_rho'][:]

# plotting
plt.close('all')
fig = plt.figure(figsize=(14,10))

vn = 'tide_Cangle'

f0 = ds0[vn][0,:,:]
f2 = ds2[vn][0,:,:]

f0[m==0] = np.nan
f2[m==0] = np.nan

ax = fig.add_subplot(131)
cs = ax.pcolormesh(f0)
fig.colorbar(cs)

ax = fig.add_subplot(132)
cs = ax.pcolormesh(f2)
fig.colorbar(cs)

ax = fig.add_subplot(133)
cs = ax.pcolormesh(f0 - f2, vmin=-0.01, vmax=0.01, cmap='bwr')
fig.colorbar(cs)

plt.show()

