#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to plot netcdf bathymetry files.
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

import os
dir0 = os.environ.get('HOME') + '/Documents/'

import os
import sys
pth = os.path.abspath(dir0 +'LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import Lfun
    
pth = os.path.abspath(dir0 +'LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

from importlib import reload
import gfun_utility as gfu
reload(gfu)

tdir0 = dir0 + 'ptools_data/topo/'

plot_dir = tdir0 + 'Figs_topo/'
Lfun.make_dir(plot_dir)

tdict = {'srtm15': ['topo15.nc'],
         'cascadia': ['cascadia_gridded.nc'],
         'psdem': ['PS_183m.nc', 'PS_91m.nc', 'PS_27m.nc'],
         'ttp_patch': ['TTP_Regional_27m_patch.nc']}

plt.close('all')

for k in tdict.keys():
    for nc_fn in tdict[k]:
        print(50*'=')
        print(nc_fn)
        out_fn = tdir0 + k + '/' + nc_fn
        # plot from NetCDF
        ds = nc.Dataset(out_fn)
        zfun.ncd(out_fn)
        lon = ds['lon'][:]
        lat = ds['lat'][:]
        z = ds['z'][:]
        ds.close()
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        plon = lon[:-1] + np.diff(lon)/2
        plat = lat[:-1] + np.diff(lat)/2
        aa = [plon[0], plon[-1], plat[0], plat[-1]]
        vv = int(np.std(z))
        cs = ax.pcolormesh(plon, plat, z[1:-1, 1:-1],
                      cmap='terrain', vmin=-vv, vmax=100)
        fig.colorbar(cs, ax=ax, extend='both')
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
        R = zfun.earth_rad(np.mean(lat))
        dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(lon[:2])/180)
        dy = R * (np.pi*np.diff(lat)/180)
        minres = np.min([np.min(dx), np.min(dy)])
        ax.set_title('%s/%s: minres = %d m' % (k, nc_fn, int(minres)))
        ax.text(.05,.05,'('+str(len(lat))+','+str(len(lon))+')', size=20, transform=ax.transAxes)
        plt.show()
    
        plt.savefig(plot_dir + nc_fn[:-3] + '.png')
