"""
Plot results of a particle tracking experiment for Dungeness Crab Sea Grant
proposal with Sean McDonald.  2021.05.11

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
sys.path.append(os.path.abspath('../plotting'))
import pfun
import seawater as sw

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

Ldir = Lfun.Lstart()

# Set an experiment to plot from.
indir0 = Ldir['LOo'] + 'tracks/'
indir = 'dcrab1_3d/'
rel = 'release_2018.03.01.nc'

# get Datasets
fn = indir0 + indir + rel
fng = indir0 + indir + 'grid.nc'
dsr = nc.Dataset(fn)
dsg = nc.Dataset(fng)

NT, NP = dsr['lon'].shape

# get a list of datetimes
ot_vec = dsr['ot'][:]
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]
hh = dsg['h'][:]
maskr = dsg['mask_rho'][:] # 1=water, 0=land
#
# tracks (packed [time, particle])
lon = dsr['lon'][:]
lat = dsr['lat'][:]
x0 = lon[0,:]
y0 = lat[0,:]

# make masks for each release location
rloc_dict = {'South Sound': (-122.7, 47.11),
            'Hood Canal': (-122.86, 47.67),
            'Whidbey': (-122.6, 48.25)}
rcolor_dict = {'South Sound': 'r',
            'Hood Canal': 'b',
            'Whidbey': 'g'}
rloc_idict = dict()
mask_arr = np.ones((3, NP)) == 0 # Bool array, initialized as all False
rcount = 0
for rloc in rloc_dict.keys():
    rloc_idict[rloc] = rcount
    xr, yr = rloc_dict[rloc]
    for ii in range(NP):
        dist, ang = sw.dist([x0[ii], xr], [y0[ii], yr])
        if dist <= 10:
            mask_arr[rcount, ii] = True
    rcount += 1

# PLOTTING
plt.close('all')
fs = 16
plt.rc('font', size=fs)
fig = plt.figure(figsize=(10,12))

# MAP
# set domain limits
aa = [-124, -122, 47, 49]
ax = fig.add_subplot(111)
# zm = -np.ma.masked_where(maskr==0, hh)
# plt.pcolormesh(lonp, latp, zm[1:-1, 1:-1], vmin=-100, vmax=0,
#     cmap='terrain', alpha=.25)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks([-124, -123, -122])
ax.set_yticks([47, 48, 49])
#ax.set_title(indir.strip('/'))

# add the tracks (packed [time, particle])
for rloc in rloc_dict.keys():
    rcount = rloc_idict[rloc]
    c = rcolor_dict[rloc]
    mask = mask_arr[rcount,:]
    Lon = lon[:,mask]
    Lat = lat[:,mask]
    ax.plot(Lon, Lat, linestyle='-',color=c, linewidth=.1, alpha=.3)
    #ax.plot(Lon[0,:], Lat[0,:], '.k', alpha=.3)
    #ax.plot(Lon[-1,:], Lat[-1,:], 'o'+c, alpha=.3)
    ax.text(rloc_dict[rloc][0], rloc_dict[rloc][1], rcount+1, color=c, weight='bold', size=fs*1.5,
        bbox=dict(facecolor='w', edgecolor=c, alpha=.6), ha='center', va='center')

#fig.tight_layout()
plt.show()
plt.rcdefaults()

dsr.close()
dsg.close()

