"""
Code to help debug gen_cmop_nudge.py
"""

import netCDF4 as nc
import numpy as np
import os, sys
from dateutil import parser
from scipy import interp
from scipy.spatial import cKDTree
from scipy.io import FortranFile
from scipy import interpolate

import xarray

from lo_tools import zfun, zrfun
from time import time

import critfc_functions as crfun

hgrid = '/Users/pm8/Documents/LO_data/critfc/hgrid.ll'
vgrid = '/Users/pm8/Documents/LO/post/critfc0/vgrid.in'
depthfile = '/Users/pm8/Documents/LiveOcean_roms/output/cas6_v3_lo8b/f2019.07.04/ocean_his_0001.nc'
ncfile = depthfile

try:
    grid = crfun.readHGrid(hgrid, True)
    x = grid.nodes[:,1]
    y = grid.nodes[:,2]
    d = grid.nodes[:,3]
    nnodes = len(x)
except:
    print('''Unable to open file %s''' % (hgrid,))
    sys.exit(1)
    
try:
    vCoords = crfun.vc.fromVGridFile(vgrid)
except Exception as e:
    print(e)
    
try:
    # construct vertical grid
    print('Vert grid')
    Z, kbp2, iwet = vCoords.computeVerticalCoordinates(d * 0.0,d)
except Exception as e:
    print(e)
    print('''Unable to open vertical grid file %s''' % (vgrid,))
    sys.exit(1)

nvrt = vCoords.nvrt

# print(hgrid, startdate, nnodes)

#print('ocean depths')
pncid = nc.Dataset(depthfile, 'r');
xrho = pncid.variables['lon_rho'][:]
yrho = pncid.variables['lat_rho'][:]
xr=xrho[:]
yr=yrho[:]

# PM Edit
#depths = pncid.variables['depths'][:]
# make the depths file
G, S, T = zrfun.get_basic_info(depthfile)
depths = -zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
depn = depths.shape[0]

# ROMS depths interpolated to SELFE grid
dep=np.zeros((depn,nnodes))
for i in range(0,depn):
    fz=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],depths[i,:,:].T,kx=1,ky=1)
    dep[i,:]=fz.ev(x,y)
    
temp_prev=np.zeros((depn,nnodes))
salt_prev=np.zeros((depn,nnodes))
print('processing')

# if testing:
#     Nfiles = 3
# else:
#     Nfiles = nfiles
    
# for i in range(Nfiles):
# there is no file #1 use #2 twice
tt0 = time()

# tfile = startf + i + 1 #add 1 because there is an offset between file # and time step
# ncfile = '''%s/ocean_his_%04d.nc''' % (basedir, tfile,)
#print(ncfile)
otime = i * 3600
#for i in range(len(x)):
temp_new=np.zeros((depn,nnodes))
salt_new=np.zeros((depn,nnodes))

print('Step %d: salt_new.max = %0.2f' % (1, np.nanmax(salt_new)))

if True:
    psncid = nc.Dataset(ncfile, 'r');
    ocean_time = psncid.variables['ocean_time'][0]
    salt = psncid.variables['salt'][0]
    temp = psncid.variables['temp'][0]
    psncid.close()
else:
    ds = xarray.open_dataset(ncfile)
    salt = ds.salt[0,:,:,:]
    temp = ds.temp[0,:,:,:]
    ds.close()


j = 29

fld = salt[j,:,:].T
ff = fld.data
ff[fld.mask] = 1e37
#ffm = np.ma.masked_where(np.isnan(ff), ff)
Xrho = xrho[0,:].data
Yrho = yrho[:,0].data
fs=interpolate.RectBivariateSpline(Xrho,Yrho,ff,kx=1,ky=1)
a=fs(x,y, grid=False)
print(np.nanmax(a))

# for j in range(0,depn):
#     ft=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],temp[j,:,:].T,kx=1,ky=1)
#     temp_new[j,:]=ft.ev(x,y)
#     fs=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],salt[j,:,:].T,kx=1,ky=1)
#     salt_new[j,:]=fs.ev(x,y)
#
# print('Step %d: salt_new.max = %0.2f' % (2, np.nanmax(salt_new)))