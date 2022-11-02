"""
Extract sections of transport, salinity, and other variables
through sections created using create_sections.py.

The new capabililty here is that the sections can be diagonal
and even made of multiple segments.

Example call:
run extract_sections -gtx cas6_v00_uu0mb -0 2021.07.04 -1 2021.07.06 -ctag test -test True

"""

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cmocean import cm
from time import time

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# set list of variables to extract
# if Ldir['get_bio']:
#     vn_list = tef_fun.vn_list
# else:
#     vn_list = ['salt']
vn_list = ['salt']

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt00 = time()
print(' Doing TEF2 extraction for '.center(60,'='))
print(' gtagex = ' + Ldir['gtagex'])
outname = 'extractions_' + ds0 + '_' + ds1
print(' outname = ' + outname)

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2' / outname
Lfun.make_dir(out_dir, clean=True)

# make the scratch directory for holding temporary files
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('tef2_temp_' + ds0 + '_' + ds1)
Lfun.make_dir(temp_dir, clean=True)

# input and output locations
grid_fn = Ldir['grid'] / 'grid.nc'
collection_name = Ldir['gridname'] + '_' + Ldir['collection_tag']
collection_dir = Ldir['LOo'] / 'extract' / 'tef2' / collection_name
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
h[m==0] = np.nan
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
Plon = plon[0,:]
Plat = plat[:,0]
aa = pfun.get_aa(ds)
ds.close

# The code below will eventually get wrapped into a tef2 version of 
# extract_sections_one_time.py

def get_sn_list():
    df_list = list(collection_dir.glob('*.p'))
    sn_list = [item.name.replace('.p','') for item in df_list]
    return sn_list

sn_list = get_sn_list()

if Ldir['testing'] == True:
    sn_list = [sn_list[0]]
    
# initialize plot
plt.close('all')
if True:
    fig = plt.figure(figsize=(8,8)) # laptop size
else:
    fig = plt.figure(figsize=(12,12)) # external monitor size
ax = fig.add_subplot(111)
ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
pfun.dar(ax)
ax.text(.05,.95,Ldir['gridname'],transform=ax.transAxes,
    fontweight='bold')
plt.show()

# dicts to keep line and text handles
ld = dict(); td = dict()

def get_xy(sn):
    df = pd.read_pickle(collection_dir / (sn + '.p'))
    x = df.x.to_numpy()
    y = df.y.to_numpy()
    return x, y
    
def get_stairstep(sn):
    # still need to deal with lines that have multiple segments
    x, y = get_xy(sn)
    NP = len(x)
    for ii in range(NP-1):
        x0 = x[ii]; y0 = y[ii]
        x1 = x[ii+1]; y1 = y[ii+1]
        ix0 = zfun.find_nearest_ind(Plon,x0)
        iy0 = zfun.find_nearest_ind(Plat,y0)
        ix1 = zfun.find_nearest_ind(Plon,x1)
        iy1 = zfun.find_nearest_ind(Plat,y1)
    ix, iy = zfun.get_stairstep(ix0, ix1, iy0, iy1)
    return ix, iy

def plot_line(sn):
    x, y = get_xy(sn)
    ax.plot(x,y,'-o',c='orange', lw=2)
    nudge = (aa[3] - aa[2]) / 200
    ax.text(x[0],y[0]-nudge,sn,c='orange',ha='center',va='top',
        fontweight='bold')
    plt.draw()
    
for sn in sn_list:
    plot_line(sn)
    ix, iy = get_stairstep(sn)
    NP = len(ix)
    for ii in range(NP-1):
        ax.plot([Plon[ix[ii]],Plon[ix[ii+1]]], [Plat[iy[ii]],Plat[iy[ii+1]]],'-or')
plt.draw()
    
    
