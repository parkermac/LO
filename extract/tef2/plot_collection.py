"""
Code to plot an existing collection of sections.

Command line arguments:

-gctag: [gridname]_[collection_tag]
-small: True to work on a laptop

Examples:

run plot_collection -gctag cas6_c0

"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cmocean import cm
import sys

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# gridname and tag for collection folder
parser.add_argument('-gctag', default='cas6_c0', type=str)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input location
gctag = args.gctag
gridname = gctag.split('_')[0]
Ldir = Lfun.Lstart(gridname=gridname)
in_dir = Ldir['LOo'] / 'extract' / 'tef2' / ('sections_' + gctag)
    
# get grid data
ds = xr.open_dataset(Ldir['grid'] / 'grid.nc')
h = ds.h.values
m = ds.mask_rho.values
h[m==0] = np.nan
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
ds.close

# initialize plot
plt.close('all')
if args.small == True:
    fig = plt.figure(figsize=(8,8)) # laptop size
else:
    fig = plt.figure(figsize=(12,12)) # external monitor size
ax = fig.add_subplot(111)
ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
pfun.dar(ax)
ax.text(.05,.95,gctag,transform=ax.transAxes,
    fontweight='bold',bbox=pfun.bbox)

def plot_line(sn):
    df = pd.read_pickle(in_dir / (sn + '.p'))
    x = df.x.to_numpy()
    y = df.y.to_numpy()
    ax.plot(x,y,'-o',c='orange', lw=2)
    ax.plot(x[0],y[0],'or',lw=2)
    ax.text(x[0],y[0],'\n'+sn,c='r',ha='center',va='top',
        fontweight='bold')
    plt.draw()

def get_sn_list():
    df_list = list(in_dir.glob('*.p'))
    sn_list = [item.name.replace('.p','') for item in df_list]
    return sn_list

# plot all existing sections in the collection
sn_list = get_sn_list()
for sn in sn_list:
    plot_line(sn)
    
plt.show()

