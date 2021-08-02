"""
Generic code to plot any mooring extraction, using pcolor.
"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import plotting_functions as pfun

#import netCDF4 as nc
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cmocean

Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
moor_name = Lfun.choose_item(in_dir, tag='.nc', exclude_tag='', itext='** Choose mooring extraction from list **')
moor_fn = in_dir / moor_name

# load everything using xarray
xs = xr.load_dataset(moor_fn)
ot = xs.ocean_time.values
ot_dt = pd.to_datetime(ot)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
T = t/86400 # time in days from start
print('time step of mooring'.center(60,'-'))
print(t[1])
print('time limits'.center(60,'-'))
print('start ' + str(ot_dt[0]))
print('end   ' + str(ot_dt[-1]))
print('info'.center(60,'-'))
VN_list = []
for vn in xs.data_vars:
    print('%s %s' % (vn, xs[vn].shape))
    VN_list.append(vn)
    
s = xs['salt'].values
th = xs['temp'].values
z = xs['z_w'].values
NT, NZ = z.shape

# make variables at midddle times, for pcolormesh
S = (s[1:,:] + s[:-1,:])/2
TH = (th[1:,:] + th[:-1,:])/2

plt.close('all')
pfun.start_plot()
fig = plt.figure()

ax = fig.add_subplot(211)
cs = ax.pcolormesh(T.reshape((NT,1))*np.ones((1,NZ)), z, S, cmap=cmocean.cm.haline)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_ylabel('Z [m]')
ax.set_title(moor_name)
ax.text(.05, .1, 'Salinity [g/kg]', c='k', weight='bold', transform=ax.transAxes)

ax = fig.add_subplot(212)
cs = ax.pcolormesh(T.reshape((NT,1))*np.ones((1,NZ)), z, TH, cmap=cmocean.cm.balance)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_xlabel('Time [days from start of record]')
ax.set_ylabel('Z [m]')
ax.text(.05, .1, 'Potential Temperature [deg C]', c='w', weight='bold', transform=ax.transAxes)

plt.show()

pfun.end_plot()