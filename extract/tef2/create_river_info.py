"""
Code to create the river info to associate with a set of segments.

This will rely on a rivers.nc file from the LO forcing.

Example command:
run create_river_info -gridname cas7 -frc trapsV00 -dstr 2017.07.04 -test True

Use -test True to see a useful plot of all the river source locations.

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
parser.add_argument('-gridname', default='cas7', type=str)
parser.add_argument('-frc', default='trapsV00', type=str)
parser.add_argument('-dstr',default='2017.07.04', type=str)
parser.add_argument('-test', '--testing', default=True, type=Lfun.boolean_string)
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

# Vertical sources, which have df.dir==2, are on the rho grid.
df.loc[(df.dir==2), 'irho'] = df.loc[(df.dir==2), 'ii']
df.loc[(df.dir==2), 'jrho'] = df.loc[(df.dir==2), 'jj']

# check results and save them if there are no bad points
ir = df.irho.to_numpy(dtype=int)
jr = df.jrho.to_numpy(dtype=int)
xr = lor[ir]
yr = lar[jr]
source_mask = m[jr,ir]
ngood = (source_mask==1).sum()
nbad = (source_mask==0).sum()
print('Number of good sources = %d' % (int(ngood)))
print('Number of bad sources = %d' % (int(nbad)))

if nbad == 0:
    print('\nSaving DataFrame to:\n%s' % (str(out_fn)))
    pd.to_pickle(df, out_fn)
else:
    print('\nNot saving results')

if testing == True:
    
    # plotting
    plt.close('all')
    pfun.start_plot()
    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(111)

    ax.pcolormesh(plon,plat,h, cmap='cool')
    pfun.dar(ax)

    # plot locations of recipient rho cells
    ax.plot(xr,yr,'ok')

    # plot source locations and directions

    # u grid
    dfup = df[(df.dir==0) & (df.sgn==1)]
    ax.plot(lou[dfup.iu.to_numpy(dtype=int)], lau[dfup.ju.to_numpy(dtype=int)],'>k')
    dfum = df[(df.dir==0) & (df.sgn==-1)]
    ax.plot(lou[dfum.iu.to_numpy(dtype=int)], lau[dfum.ju.to_numpy(dtype=int)],'<k')

    # v grid
    dfvp = df[(df.dir==1) & (df.sgn==1)]
    ax.plot(lov[dfvp.iv.to_numpy(dtype=int)], lav[dfvp.jv.to_numpy(dtype=int)],'^k')
    dfvm = df[(df.dir==1) & (df.sgn==-1)]
    ax.plot(lov[dfvm.iv.to_numpy(dtype=int)], lav[dfvm.jv.to_numpy(dtype=int)],'vk')

    # rho grid
    dfrp = df[(df.dir==2) & (df.sgn==1)]
    ax.plot(lor[dfrp.irho.to_numpy(dtype=int)], lar[dfrp.jrho.to_numpy(dtype=int)],'or')
    dfrm = df[(df.dir==2) & (df.sgn==-1)] # -1 should not exist for vertical sources
    ax.plot(lor[dfrm.irho.to_numpy(dtype=int)], lar[dfrm.jrho.to_numpy(dtype=int)],'og')
    
    # add names
    for II in df.index:
        ax.text(lor[int(df.loc[II,'ii'])],
                lar[int(df.loc[II,'jj'])],
                df.loc[II,'name'])

    plt.show()

