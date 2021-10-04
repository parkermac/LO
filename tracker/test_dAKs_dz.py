"""
Code to test calculation of dAKs/dz for vmix.
"""

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import zfun

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# get fields
fn = ('/Users/pm8/Documents/LiveOcean_output/moor/' +
    'cas6_v3_lo8b_2019.07.04_2019.07.04/vmixWB_hourly.nc')
ds = nc.Dataset(fn)
ii = 24
k = ds['AKs'][ii,:].data
zw = ds['z_w'][ii,:].data
zr = ds['z_rho'][ii,:].data
ds.close()

# adjust z's to go to zero at the free surface
ssh = zw[-1]
zw = zw - ssh
zr = zr - ssh
z0 = zw[0] # z at the bottom

# calculations
dz = np.diff(zw)
dk = np.diff(k)
dkdz = dk/dz # on the rho-grid

if True:
    # smooth k a bit
    kf = zfun.filt_hanning(k, n=3)
    kf[0] = k[0]
    kf[-1] = k[-1]
    k = kf
    dk = np.diff(k)
    dkdz = dk/dz # on the rho-grid

def get_kp_dkdzp(zp, zw, zr, k, dkdz):
    do_nearest = True
    # find k at particle positions
    i0, i1, fr = zfun.get_interpolant(zp, zw, extrap_nan=False)
    if do_nearest:
        fr = np.round(fr)
    kp = (1-fr)*k[i0] + fr*k[i1]
    # find dkdz at particle positions
    i0, i1, fr = zfun.get_interpolant(zp, zr, extrap_nan=False)
    if do_nearest:
        fr = np.round(fr)
    dkdzp = (1-fr)*dkdz[i0] + fr*dkdz[i1]
    return kp, dkdzp

# initialize particles
NP = 4000
zp = np.linspace(z0,0,NP)

# move them one or more time steps
dt = 300 # seconds
ndays = 1
NT = int(ndays*86400/dt)
for tt in range(NT): # run for a ndays
    
    # first move partway to get appropriate k and dkdz
    kp, dkdzp = get_kp_dkdzp(zp, zw, zr, k, dkdz)
    dz_part = dt*dkdzp/2
    kp, dkdzp = get_kp_dkdzp(zp + dz_part, zw, zr, k, dkdz)
    
    # then diffuse the particles for real
    rand = np.random.standard_normal(NP)
    w = rand*np.sqrt(2*kp/dt) + dkdzp
    zp += w*dt
    
    # enforce limits using reflection
    mask = zp > 0
    zp[mask] = -zp[mask]
    mask = zp < z0
    zp[mask] = z0 + (z0-zp[mask])

# plotting
plt.close('all')
fig = plt.figure(figsize=(12,12))

ax = fig.add_subplot(131)
ax.plot(k, zw, '-ob', lw=2)
ax.set_title('AKs')
ax.set_ylim(z0,0)

ax = fig.add_subplot(132)
ax.plot(dkdz, zr, '-or', lw=2)
ax.set_title('dAKs/dz')
ax.set_ylim(z0,0)

ax = fig.add_subplot(133)
if False:
    ax.plot(np.random.random(NP), zp, 'og', ms=2)
    ax.set_title('Particles')
    ax.set_xlim(0,1)
    ax.set_ylim(z0,0)
else:
    NB = 28 # number of bins (assumes NP = 4000)
    amin = 102.5 * np.ones(2)
    amax = 187.3 * np.ones(2)
    bins=np.linspace(z0, 0, NB+1)
    cbins = bins[:-1] + np.diff(bins)/2
    counts, obins = np.histogram(zp, bins=bins)
    ax.plot(counts, cbins, '-og', lw=2)
    ax.plot(amin, [z0,0],'-k',lw=3, alpha=.3)
    ax.plot(amax, [z0,0],'-k',lw=3, alpha=.3)
    ax.set_ylim(z0,0)
    ax.set_xlim(0, 300)
    ax.set_title('Counts')

plt.show()