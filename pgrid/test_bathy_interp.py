"""
Code to test bathymetry interpolation and overlay.
"""

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

from pathlib import Path
import xarray as xr
import netCDF4 as nc
import gfun_utility as gfu
import numpy as np
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart()

aa = [-124, -122, 47, 49]
res = 600 # target resolution (m)
Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
lon, lat = np.meshgrid(Lon_vec, Lat_vec)
plon, plat = pfun.get_plon(lon, lat)

t_dir = Ldir['data'] / 'topo'
# list of topo files: coarsest to finest
t_list = [t_dir / 'srtm15' / 'topo15.nc',
          t_dir / 'cascadia' / 'cascadia_gridded.nc',
         t_dir / 'psdem' / 'PS_183m.nc',
         t_dir / 'ttp_patch' / 'TTP_Regional_27m_patch.nc']
#t_list = [t_dir / 'psdem' / 'PS_183m.nc']

plt.close('all')
for t in t_list:

    fn = t_dir / t
    if False:
        ds = nc.Dataset(fn)
        tlon_vec = ds['lon'][:]
        tlat_vec = ds['lat'][:]
        UU = ds['z'][:]
        U = zfun.fillit(UU)
    else:
        ds = xr.open_dataset(fn)
        tlon_vec = ds.lon.values
        tlat_vec = ds.lat.values
        U = ds.z.values
        # There is a bug in xarray with these files: it does
        # not set masked regions to nan.  So we do it by hand.
        U[U>1e6] = np.nan
    ds.close()
    
    if False:
        xi0, xi1, xf = zfun.get_interpolant(Lon_vec, tlon_vec, extrap_nan=True)
        yi0, yi1, yf = zfun.get_interpolant(Lat_vec, tlat_vec, extrap_nan=True)
        NR = len(Lat_vec)
        NC = len(Lon_vec)
        XF = xf.reshape((1,NC)) * np.ones((NR,1))
        YF = yf.reshape((NR,1)) * np.ones((1,NC))
        # bi linear interpolation
        u00 = U[yi0,:][:,xi0]
        u10 = U[yi1,:][:,xi0]
        u01 = U[yi0,:][:,xi1]
        u11 = U[yi1,:][:,xi1]
        ui = (1-YF)*((1-XF)*u00 + XF*u01) + YF*((1-XF)*u10 + XF*u11)
    else:
        tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
        ui = zfun.interp2(lon, lat, tlon, tlat, U)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    cs = ax.pcolormesh(plon, plat, ui)

    plt.show()


