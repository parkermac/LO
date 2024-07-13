"""
Plot results of a particle tracking experiment.

Development version.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

# Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

# get Datasets
fn = in_dir0 / exp_name / rel
fng = in_dir0 / exp_name / 'grid.nc'
dsr = xr.open_dataset(fn, decode_times=False)
# dsg = xr.open_dataset(fng)

NT, NP = dsr.lon.shape

# get a list of datetimes
ot_vec = dsr.ot.values
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# # gather some fields, for convenience
# lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
# hh = dsg.h.values
# maskr = dsg.mask_rho.values
#

# # subsample output for plotting
# npmax = 300 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))

# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]

# # make a mask that is False from the time a particle first leaves the domain
# # and onwards
# AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
#         dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
# ib_mask = np.ones(lon.shape, dtype=bool)
# ib_mask[lon < AA[0]] = False
# ib_mask[lon > AA[1]] = False
# ib_mask[lat < AA[2]] = False
# ib_mask[lat > AA[3]] = False
# NTS, NPS = lon.shape
# for pp in range(NPS):
#     tt = np.argwhere(ib_mask[:,pp]==False)
#     if len(tt) > 0:
#         ib_mask[tt[0][0]:, pp] = False

# # and apply the mask to lon and lat
# lon[~ib_mask] = np.nan
# lat[~ib_mask] = np.nan

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(14,8))

# # MAP
# # set domain limits
# if False:
#     # plot full domain
#     aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
# else:
#     # automatically plot region of particles, with padding
#     pad = .02
#     aa = [np.nanmin(lon) - pad, np.nanmax(lon) + pad,
#     np.nanmin(lat) - pad, np.nanmax(lat) + pad]
    
# ax = fig.add_subplot(121)
# zm = -np.ma.masked_where(maskr==0, hh)
# plt.pcolormesh(lonp, latp, zm, vmin=-100, vmax=0,
#     cmap='terrain', alpha=.25)
# pfun.add_coast(ax)
# ax.axis(aa)
# pfun.dar(ax)
# ax.set_xlabel('Longitude')
# ax.set_ylabel('Latitude')
# ax.set_title(exp_name.strip('/'))
# # add the tracks (packed [time, particle])
# # regular spaghetti plots
# ax.plot(lon, lat, '-k', linewidth=.2)
# ax.plot(lon[0,:], lat[0,:], 'og', alpha=.3)
# ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)

for ii in range(NP):

    fig = plt.figure()
    # time series
    ts = ot_vec - ot_vec[0]
    td = ts/86400
    h = dsr.h[:,ii].values
    z = dsr.z[:,ii].values
    cs = dsr.cs[:,ii].values
    w = dsr.w[:,ii].values
    lat = dsr.lat[:,ii].values
    zeta = dsr.zeta[:,ii].values

    NR = 2

    ax = fig.add_subplot(NR,1,1)
    ax.plot(td,z,'-r')
    ax.plot(td,-h,'-k')
    ax.plot(td,zeta,'-b')
    ax.set_ylabel('z')
    ax.set_title('Particle number %d' % (ii))

    ax = fig.add_subplot(NR,1,2)
    ax.plot(td,cs,'-r')
    ax.set_ylabel('cs')


plt.show()
pfun.end_plot()

dsr.close()
# dsg.close()

