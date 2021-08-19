"""
Code to compare the results of atm0 to the LoveOcean atm1, to make sure
it is all working in LO.

RESULT: The new version results exactly match the old ones.
"""

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

date_string = '2019.07.04'

Ldir = Lfun.Lstart(gridname='cas6', tag='v3')

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle

import atm_fun as afun

# get files

grid_fn = str(Ldir['grid']) + '/grid.nc'
dsg = nc.Dataset(grid_fn)
Lon = dsg['lon_rho'][:]
Lat = dsg['lat_rho'][:]
dsg.close()

# new
indir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + date_string) / 'atm0'

#old
indiro = Ldir['parent'] / 'LiveOcean_output' / 'cas6_v3' / ('f' + date_string) / 'atm1'
    
#%% plotting
aa = [-130, -122, 42, 52]
lim_dict = dict(zip(afun.outvar_list, afun.lim_list))

plt.close('all')

if False:
    outvar_list = ['Uwind'] # testing
else:
    outvar_list = afun.outvar_list

for vn in outvar_list:
    
    ds = nc.Dataset(indir / (vn + '.nc'))
    dso = nc.Dataset(indiro / (vn + '.nc'))
    
    vmin = lim_dict[vn][0]
    vmax = lim_dict[vn][1]

    # extract variables to plot
    forecast_hour = 20
    vv = ds[vn][forecast_hour, :, :]
    vvo = dso[vn][forecast_hour, :, :]
    
    ds.close()
    dso.close()
        
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,7), squeeze=False)

    ax = axes[0,0]
    cs = ax.pcolormesh(Lon, Lat, vvo, cmap='rainbow', vmin=vmin, vmax=vmax, shading='nearest')
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('Hour = ' + str(forecast_hour) + ', LO atm0')
    ax.text(.97, .15, vn,
        horizontalalignment='right', transform=ax.transAxes,
        fontsize=18, fontweight='bold')
    ax.text(.97, .10, Ldir['gtag'],
        horizontalalignment='right', transform=ax.transAxes)
    ax.text(.97, .05, date_string,
        horizontalalignment='right', transform=ax.transAxes)
    ax.axis(aa)

    ax = axes[0,1]
    cs = ax.pcolormesh(Lon, Lat, vv, cmap='rainbow', vmin=vmin, vmax=vmax, shading='nearest')
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('LiveOcean atm1')
    ax.axis(aa)

    ax = axes[0,2]
    dvv = vv - vvo
    dv = np.max(np.abs([np.max(dvv), np.min(dvv)]))
    cs = ax.pcolormesh(Lon, Lat, dvv, cmap='bwr', vmin=-dv, vmax=dv, shading='nearest')
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('New - Old')
    ax.axis(aa)
    
    fig.tight_layout()

plt.show()

