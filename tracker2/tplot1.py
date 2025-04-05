"""
Plot results of a particle tracking experiment.

Like tplot.py but plots all particles and no tracks.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

# Choose an experiment and release to plot.
in_dir00 = Ldir['LOo'] / 'tracks2'
gtx_name = Lfun.choose_item(in_dir00, tag='', exclude_tag='.csv',
    itext='** Choose gtagex from list **', last=False)
in_dir0 = in_dir00 / gtx_name
exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

# get Datasets
fn = in_dir0 / exp_name / rel
fng = in_dir0 / exp_name / 'grid.nc'
dsr = xr.open_dataset(fn, decode_times=False)
dsg = xr.open_dataset(fng)

lon0 = dsr.lon[0,:].to_numpy()
lat0 = dsr.lat[0,:].to_numpy()
lon1 = dsr.lon[-1,:].to_numpy()
lat1 = dsr.lat[-1,:].to_numpy()

NT, NP = dsr.lon.shape

# get a list of datetimes
ot_vec = dsr.ot.values
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
hh = dsg.h.values
maskr = dsg.mask_rho.values

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(10,12))

fig = plt.figure()

# MAP
# set domain limits
# automatically plot region of particles, with padding
pad = .5
aa = [np.nanmin(lon0) - pad, np.nanmax(lon0) + pad,
np.nanmin(lat0) - pad, np.nanmax(lat0) + pad]
    
ax = fig.add_subplot(111)
zm = -np.ma.masked_where(maskr==0, hh)
plt.pcolormesh(lonp, latp, zm, vmin=-100, vmax=0,
    cmap='terrain', alpha=.25)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(exp_name.strip('/'))
ax.plot(lon0, lat0, '.g', alpha=.3)
ax.plot(lon1, lat1, '.r', alpha=.2)

plt.show()
pfun.end_plot()

dsr.close()
dsg.close()

