"""
Code to create the river info to associate with a set of segments.

This will rely on a rivers.nc file from the LO forcing.
"""

from lo_tools import Lfun, zrfun, zfun
from lo_tools import plotting_functions as pfun
from time import time
import sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cmocean import cm

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# info to find a rivers.nc file
parser.add_argument('-gridname', default='cas6', type=str)
parser.add_argument('-frc', default='traps2', type=str)
parser.add_argument('-dstr',default='2017.07.04', type=str)

parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

gridname = args.gridname
dstr = args.dstr
frc = args.frc

# input and output locations
testing = args.testing
Ldir = Lfun.Lstart(gridname=gridname)

grid_fn = Ldir['grid'] / 'grid.nc'
riv_fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + dstr) / frc / 'rivers.nc'

out_name = 'riv_df_' + gridname + '_' + frc + '.p'
out_dir = Ldir['LOo'] / 'extract' / 'tef2'
out_fn = out_dir / out_name
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
# depth for plotting
h[m==0] = np.nan
# coordinates for plotting
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
# coordinates for convenience in plotting
lor = ds.lon_rho[0,:].values
lar = ds.lat_rho[:,0].values
lou = ds.lon_u[0,:].values
lau = ds.lat_u[:,0].values
lov = ds.lon_v[0,:].values
lav = ds.lat_v[:,0].values
ds.close

# get river info
rds = xr.open_dataset(riv_fn)

# Pull out river info into a pandas DataFrame
df = pd.DataFrame(columns=['name','ii','jj','dir','sgn', 'iu','ju', 'iv','jv', 'irho','jrho'])
df.loc[:,'name'] = list(rds.river_name.values)
df.loc[:,'ii'] = list(rds.river_Xposition.values)
df.loc[:,'jj'] = list(rds.river_Eposition.values)
df.loc[:,'dir'] = list(rds.river_direction.values) # 0 = u-grid, 1 = v-grid, 2 = rho-grid

# also need to know if the source is positive or negative in order to assign it to the
# rho grid point it flows into
q = rds.river_transport.values
df.loc[:,'sgn'] = np.sign(q.mean(axis=0))

# translate u-grid sources to dataset indices
df.loc[df.dir==0,'iu'] = df.loc[df.dir==0,'ii'] - 1
df.loc[df.dir==0,'ju'] = df.loc[df.dir==0,'jj']

# translate v-grid sources to dataset indices
df.loc[df.dir==1,'iv'] = df.loc[df.dir==1,'ii']
df.loc[df.dir==1,'jv'] = df.loc[df.dir==1,'jj'] - 1

# use the sign of the source to figure out which rho grid point to associate with the source
# i.e. which rho point does it flow into

df.loc[(df.dir==0) & (df.sgn==1), 'irho'] = df.loc[(df.dir==0) & (df.sgn==1), 'iu'] + 1
df.loc[(df.dir==0) & (df.sgn==1), 'jrho'] = df.loc[(df.dir==0) & (df.sgn==1), 'ju'] 

df.loc[(df.dir==0) & (df.sgn==-1), 'irho'] = df.loc[(df.dir==0) & (df.sgn==-1), 'iu']
df.loc[(df.dir==0) & (df.sgn==-1), 'jrho'] = df.loc[(df.dir==0) & (df.sgn==-1), 'ju']

df.loc[(df.dir==1) & (df.sgn==1), 'irho'] = df.loc[(df.dir==1) & (df.sgn==1), 'iv']
df.loc[(df.dir==1) & (df.sgn==1), 'jrho'] = df.loc[(df.dir==1) & (df.sgn==1), 'jv'] + 1

df.loc[(df.dir==1) & (df.sgn==-1), 'irho'] = df.loc[(df.dir==1) & (df.sgn==-1), 'iv']
df.loc[(df.dir==1) & (df.sgn==-1), 'jrho'] = df.loc[(df.dir==1) & (df.sgn==-1), 'jv']

# plotting
plt.close('all')
pfun.start_plot()
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)

ax.pcolormesh(plon,plat,h, cmap='cool')
pfun.dar(ax)

ax.plot(lor[df.irho.to_numpy(dtype=int)], lar[df.jrho.to_numpy(dtype=int)],'ok')

dfu = df[df.dir==0]
ax.plot(lou[dfu.iu.to_numpy(dtype=int)], lar[dfu.ju.to_numpy(dtype=int)],'+k')

dfv = df[df.dir==1]
ax.plot(lov[dfv.iv.to_numpy(dtype=int)], lav[dfv.jv.to_numpy(dtype=int)],'+k')

plt.show()

