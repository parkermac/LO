"""
This makes a nice plot of a mooring extraction from a location in the landward
end of Admiralty Inlet, mniminally the site of the AI8 mooring in Geyer and Cannon (1982).

I put this in the tef/pm_analysis directory instead of moor because it wants to be compared
with the results of section ai4 as viewed in Qprism_series.py.
"""

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()
moor_fn = Ldir['LOo'] / 'extract' / 'cas6_v3_lo8b' / 'moor' / 'AI8_2018.01.01_2018.12.31.nc'

ds = xr.open_dataset(moor_fn)

# for vn in ds.data_vars:
#     print('%s: %s' % (vn, (str(ds[vn].shape))))
    
def rot_vec(u,v,theta):
    ur = u*np.cos(theta) + v*np.sin(theta)
    vr = v*np.cos(theta) - u*np.sin(theta)
    return ur, vr
    
def de_mean(u, v):
    up = u - np.nanmean(u)
    vp = v - np.nanmean(v)
    return up, vp
    
u = ds.u.values
v = ds.v.values
s = ds.salt.values
z_w = ds.z_w.values
dz = np.diff(z_w, axis=1)
ubar = np.sum(u*dz, axis=1) / (ds.zeta.values + ds.h.values)
vbar = np.sum(v*dz, axis=1) / (ds.zeta.values + ds.h.values)
sbar = np.sum(s*dz, axis=1) / (ds.zeta.values + ds.h.values)

# make low-passed sprime
nt, nz = s.shape
sp = s - sbar.reshape((nt,1))*np.ones((1,nz))
sp_lp = zfun.lowpass(sp, f='godin')

# find angle of principal axes
up, vp = de_mean(ubar, vbar)
theta = 0.5 * np.arctan2(2*np.nanmean(up*vp),(np.nanvar(up)-np.nanvar(vp)))

# and rotate
ubar_r, vbar_r = rot_vec(ubar, vbar, theta)
u_r, v_r = rot_vec(u,v,theta)

# and low-pass
uu = zfun.lowpass(u_r, f='godin')

# mooring location
x = ds.lon_rho.values
y = ds.lat_rho.values

# time
ot = ds.ocean_time.values
ot_dt = pd.to_datetime(ot)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
T = t/86400 # time in days from start

NT, NZ = z_w.shape
Z = z_w.mean(axis=0)

# coordinate arrays for plotting
TT = T.reshape((NT,1))*np.ones((1,NZ))
ZZ = Z.reshape((1,NZ))*np.ones((NT,1))

# make variables at middle times, for pcolormesh
U = (uu[1:,:] + uu[:-1,:])/2
S = (sp_lp[1:,:] + sp_lp[:-1,:])/2

# tidally-averaged rms velocity
Urms = np.sqrt(zfun.lowpass(ubar_r**2, f='godin'))
    
ds.close()

# PLOTTING
plt.close('all')

# 1 Check on coordinate rotation and mooring location
pfun.start_plot(figsize=(16,8))
fig = plt.figure()

ax = fig.add_subplot(121)
ax.plot(ubar, vbar, '.g', alpha=.3)
ax.plot(ubar_r, vbar_r, '.b')
ax.grid(True)
ax.axis('square')
ax.axhline()
ax.axvline()
ax.set_xlabel('Ubar (m/s)')
ax.set_ylabel('Vbar (m/s)')
ax.text(.05, .9, 'Original', c='g', weight='bold', transform=ax.transAxes)
ax.text(.05, .8, 'Rotated', c='b', weight='bold', transform=ax.transAxes)

ax = fig.add_subplot(122)
pfun.add_coast(ax)
ax.plot(x, y, '*r')
pad = .3
ax.axis([x-pad, x+pad, y-pad, y+pad])
pfun.dar(ax)
ax.set_xlabel('Longitude (deg)')
ax.set_ylabel('Latitude (deg)')
ax.text(.05, .9, 'AI8 Mooring Location', c='r', weight='bold', transform=ax.transAxes)

plt.show()
pfun.end_plot()

# 2 Time series
pfun.start_plot(figsize=(16,8))
fig = plt.figure()

ax = fig.add_subplot(311)
ax.plot(T, Urms, '-b', lw=2)
ax.set_xlim(0,365)
ax.grid(True)
ax.set_xticklabels([])
ax.set_ylabel('RMS Ubar (m/s)')

ax = fig.add_subplot(312)
cs = ax.pcolormesh(TT, ZZ, U, cmap='RdYlBu_r', vmin=-.3, vmax=.3)
ax.set_xlim(0,365)
ax.grid(True)
ax.set_xticklabels([])
ax.set_ylabel(r'$U_{lowpass}\ (m/s)$')
ax.text(.95,.5,'Range = -0.3 to 0.3', ha = 'right', weight='bold', transform=ax.transAxes)

ax = fig.add_subplot(313)
cs = ax.pcolormesh(TT, ZZ, S, cmap='RdYlBu_r', vmin=-.7, vmax=.7)
ax.set_xlim(0,365)
ax.grid(True)
ax.set_ylabel(r'$S^{\prime}_{lowpass}\ (g/kg)$')
ax.set_xlabel('Yearday 2018')
ax.text(.95,.5,'Range = -0.7 to 0.7', ha = 'right', weight='bold', transform=ax.transAxes)

plt.show()
pfun.end_plot()
