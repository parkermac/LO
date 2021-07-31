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

import netCDF4 as nc
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

ds = nc.Dataset(moor_fn)

ot = ds['ocean_time'][:]
print('time step of mooring'.center(60,'-'))
print(np.diff(ot))

tind = [Lfun.modtime_to_datetime(tt) for tt in ot]
print('time limits'.center(60,'-'))
print('start ' + str(tind[0]))
print('end   ' + str(tind[-1]))

print('info'.center(60,'-'))
VN_list = []
for vn in ds.variables:
    print('%s %s' % (vn, ds[vn].shape))
    VN_list.append(vn)
    
t = (ot - ot[0])/86400 # time in days
s = ds['salt'][:].data
th = ds['temp'][:].data
z = ds['z_w'][:].data
NT, NZ = z.shape

# make variables at midddle times, fo pcolormesh
S = (s[1:,:] + s[:-1,:])/2
TH = (th[1:,:] + th[:-1,:])/2

plt.close('all')
pfun.start_plot()
fig = plt.figure()

ax = fig.add_subplot(211)
cs = ax.pcolormesh(t.reshape((NT,1))*np.ones((1,NZ)), z, S, cmap=cmocean.cm.haline)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_ylabel('Z [m]')
ax.set_title(moor_name)
ax.text(.05, .1, 'Salinity [g/kg]', c='k', weight='bold', transform=ax.transAxes)

ax = fig.add_subplot(212)
cs = ax.pcolormesh(t.reshape((NT,1))*np.ones((1,NZ)), z, TH, cmap=cmocean.cm.balance)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_xlabel('Time [days from start of record]')
ax.set_ylabel('Z [m]')
ax.text(.05, .1, 'Potential Temperature [deg C]', c='w', weight='bold', transform=ax.transAxes)

plt.show()

pfun.end_plot()