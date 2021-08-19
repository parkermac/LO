"""
Code to plot the various WRF grids, with a selected ROMS grid.

Assumes we are running on my mac.
"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart(gridname='cas6')

import netCDF4 as nc
import matplotlib.pyplot as plt
import seawater as sw

# where are files located, and other situational choices
wrf_dir = str(Ldir['data']) + '/wrf/'
Ldir['date_string'] = '2019.07.04'

# create list of input files
d_str = Ldir['date_string'].replace('.','')
in_dir = wrf_dir + d_str + '00/'
hr = 0 # grid info is only in the zero hour
hr_str = ('0' + str(hr))[-2:]
d2_fn = (in_dir + 'wrfout.ocean_d2.' + d_str + '00.f' + hr_str + '.0000')
d3_fn = (in_dir + 'wrfout.ocean_d3.' + d_str + '00.f' + hr_str + '.0000')
d4_fn = (in_dir + 'wrfout.ocean_d4.' + d_str + '00.f' + hr_str + '.0000')

# get model grid
gds = nc.Dataset(str(Ldir['grid']) + '/grid.nc')
lon = gds['lon_rho'][:]
lat = gds['lat_rho'][:]
gds.close()

# get WRF grid(s)
def get_wrf_grid(fn):
    wds = nc.Dataset(fn)
    lon = wds['XLONG'][:].squeeze()
    lat = wds['XLAT'][:].squeeze()
    if False:
        print('\n** ' + fn.split('/')[-1])
        vn_list = []
        for vn in wds.variables:
            vn_list.append(vn)
        print(vn_list)
    wds.close()
    # grid size info
    NR, NC = lon.shape
    jj = int(NR/2); ii = int(NC/2)
    dx_km, dd_deg = sw.dist(lat[jj,ii], [lon[jj,ii], lon[jj+1,ii+1]])
    return lon, lat, dx_km
# lat, lon are only in first file of the day
lon2, lat2, dx2_km = get_wrf_grid(d2_fn)
lon3, lat3, dx3_km = get_wrf_grid(d3_fn)
lon4, lat4, dx4_km = get_wrf_grid(d4_fn)

# PLOTTING
plt.close('all')
fs=16
plt.rc('font', size=fs)

fig = plt.figure(figsize=(16,10))

# nice plot of grid domains and spacing
def draw_box(ax, x,y, s='-', c='k', w=3, a=1, t=''):
    ax.plot(x[:,0],y[:,0], linestyle=s, color=c, linewidth=w, alpha=a)
    ax.plot(x[:,-1],y[:,-1], linestyle=s, color=c, linewidth=w, alpha=a)
    ax.plot(x[0,:],y[0,:], linestyle=s, color=c, linewidth=w, alpha=a)
    ax.plot(x[-1,:],y[-1,:], linestyle=s, color=c, linewidth=w, alpha=a)
    ax.text(x[-1,-1]-1,y[-1,-1]-.5, t, fontweight='bold', color=c,
        horizontalalignment='right', verticalalignment='top')

ax = plt.subplot2grid((1,3), (0,0), colspan=2)
draw_box(ax, lon, lat, c='b', t='LiveOcean')
draw_box(ax, lon2, lat2, c='g', t=('WRF d2\n(%0.1f km)' % (dx2_km)))
draw_box(ax, lon3, lat3, c='orange', t=('WRF d3\n(%0.1f km)' % (dx3_km)))
draw_box(ax, lon4, lat4, c='r', t=('WRF d4\n(%0.1f km)' % (dx4_km)))
ax.set_title('(a) Grid extents')
pfun.dar(ax)
pfun.add_coast(ax, color='gray')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_yticks([40, 45, 50])

ax = plt.subplot2grid((1,3), (0,2), colspan=1)
al=.5
ax.plot(lon, lat, 'ob', markersize=2, alpha=al)
ax.plot(lon2, lat2, 'sg', markersize=20, alpha=al)
ax.plot(lon3, lat3, linestyle='', marker='o', color='orange', markersize=10, alpha=al)
ax.plot(lon4, lat4, '*r', markersize=5, alpha=al)
ax.axis([-123, -122.7, 47.55, 47.9])
ax.set_title('(b) Close-up of grid spacing')
pfun.dar(ax)
pfun.add_coast(ax, color='gray', linewidth=3)
ax.set_xlabel('Longitude')
ax.set_xticks([-122.9, -122.7])
ax.set_yticks([47.6, 47.7, 47.8, 47.9])
fig.tight_layout()
plt.show()
plt.rcdefaults()