"""
Code to compare the results of ocn03 to ocn02.

RESULT: The fields look reasonabley similar except for the offest
in zeta, which was the problem we were trying to solve. So all looks
good.

"""

import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

dt0 = datetime(2024,11,23)
dstr0 = dt0.strftime(Lfun.ds_fmt)

in_dir0 = Ldir['LOo'] / 'forcing' / 'cas7' / ('f' + dstr0)
fn2 = in_dir0 / 'ocn02' / 'ocean_clm.nc'
fn3 = in_dir0 / 'ocn03' / 'ocean_clm.nc'
ds2 = xr.open_dataset(fn2)
ds3 = xr.open_dataset(fn3)

g_fn = Ldir['data'] / 'grids' / 'cas7' / 'grid.nc'
dsg = xr.open_dataset(g_fn)
plon, plat = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)

plt.close('all')
pfun.start_plot(figsize=(24,10))

vn_list = ['zeta', 'salt', 'temp', 'u', 'v', 'ubar', 'vbar']

# vn_list = ['NO3',
#  'NH4',
#  'chlorophyll',
#  'phytoplankton',
#  'zooplankton',
#  'LdetritusN',
#  'SdetritusN',
#  'LdetritusC',
#  'SdetritusC',
#  'TIC',
#  'alkalinity',
#  'oxygen']

for vn in vn_list:

    fig = plt.figure()

    slev = -1 # 0 = bottom, -1 = surface
    tlev = -1 # 0 = first time, -1 = last time
    if vn in ['zeta','ubar','vbar']:
        fld2 = ds2[vn][tlev,:,:].values
        fld3 = ds3[vn][tlev,:,:].values
    else:
        fld2 = ds2[vn][tlev,slev,:,:].values
        fld3 = ds3[vn][tlev,slev,:,:].values

    if vn in ['u', 'ubar']:
        fr2 = np.nan * np.ones((dsg.lon_rho.values).shape)
        fr2[:,1:-1] = (fld2[:,1:] + fld2[:,:-1])/2
        fld2 = fr2
        fr3 = np.nan * np.ones((dsg.lon_rho.values).shape)
        fr3[:,1:-1] = (fld3[:,1:] + fld3[:,:-1])/2
        fld3 = fr3

    if vn in ['v', 'vbar']:
        fr2 = np.nan * np.ones((dsg.lon_rho.values).shape)
        fr2[1:-1,:] = (fld2[1:,:] + fld2[:-1,:])/2
        fld2 = fr2
        fr3 = np.nan * np.ones((dsg.lon_rho.values).shape)
        fr3[1:-1,:] = (fld3[1:,:] + fld3[:-1,:])/2
        fld3 = fr3
        
    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(plon, plat, fld2)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    ax.set_title('ocn02')
    ax.text(.05,.9,vn,transform=ax.transAxes, fontweight='bold', bbox=pfun.bbox)

    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(plon, plat, fld3)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    ax.set_title('ocn03')

    ax = fig.add_subplot(133)
    cs = ax.pcolormesh(plon, plat, fld3-fld2)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    ax.set_title('03 - 02')

plt.show()


