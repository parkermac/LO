"""
Extract sections of transport, salinity, and other variables
through sections created using create_sections.py.

The new capabililty here is that the sections can be diagonal
and even made of multiple segments.

Example call:
run extract_sections -gtx cas6_v00_uu0mb -0 2021.07.04 -1 2021.07.06 -ctag c0 -test True

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
# temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('tef2_temp_' + ds0 + '_' + ds1)
# Lfun.make_dir(temp_dir, clean=True)
# Do we use this for anything?

# input and output locations
grid_fn = Ldir['grid'] / 'grid.nc'
collection_name = Ldir['gridname'] + '_' + Ldir['collection_tag']
collection_dir = Ldir['LOo'] / 'extract' / 'tef2' / collection_name
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
h[m==0] = np.nan
# coordinates for plotting
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
# coordinates for stairstep generation
lop = ds.lon_psi[0,:].values
lap = ds.lat_psi[:,0].values
lor = ds.lon_rho[0,:].values
lar = ds.lat_rho[:,0].values
lou = ds.lon_u[0,:].values
lau = ds.lat_u[:,0].values
lov = ds.lon_v[0,:].values
lav = ds.lat_v[:,0].values
ds.close

# The code below will eventually get wrapped into a tef2 version of 
# extract_sections_one_time.py

def get_sn_list():
    df_list = list(collection_dir.glob('*.p'))
    sn_list = [item.name.replace('.p','') for item in df_list]
    return sn_list

sn_list = get_sn_list()

if Ldir['testing'] == True:
    sn_list = ['mb8','mb9']
    
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
    
def get_dist(xy):
    d = []
    for ii in range(len(xy)-1):
        d.append(np.abs(xy[ii][0]-xy[ii+1][0]) + np.abs(xy[ii][1]-xy[ii+1][1]))
    d = np.array(d)
    # add a fake last distance to make it the same distance as xy
    d = np.concatenate((d,np.array([1]))).astype(float)
    return d
    
def get_stairstep(sn):
    x, y = get_xy(sn)
    NP = len(x)
    ix = np.array([])
    iy = np.array([])
    for ii in range(NP-1):
        x0 = x[ii]; y0 = y[ii]
        x1 = x[ii+1]; y1 = y[ii+1]
        ix0 = zfun.find_nearest_ind(lop,x0)
        iy0 = zfun.find_nearest_ind(lap,y0)
        ix1 = zfun.find_nearest_ind(lop,x1)
        iy1 = zfun.find_nearest_ind(lap,y1)
        this_ix, this_iy = zfun.get_stairstep(ix0, ix1, iy0, iy1)
        ix = np.concatenate((ix,this_ix))
        iy = np.concatenate((iy,this_iy))
    ix = ix.astype(int)
    iy = iy.astype(int)
    return ix, iy
    
def trim(ix, iy):
    # keep only unique entries
    xy = list(zip(ix,iy))
    seen = set()
    unique_xy = []
    for item in xy:
        if item not in seen:
            unique_xy.append(item)
            seen.add(item)
    xy = np.array(unique_xy)
    # and then clip stray ends
    d = get_dist(xy)
    mask = d > 1
    trimmed = 0
    while sum(mask) >= 1:
        xy = xy[~mask]
        d = get_dist(xy)
        mask = d > 1
        trimmed += 1
    ix = [item[0] for item in xy]
    iy = [item[1] for item in xy]
    return ix, iy, trimmed

def plot_line(sn):
    x, y = get_xy(sn)
    ax.plot(x,y,'-o',c='orange', lw=2)
    ax.text(x[0],y[0],'\n'+sn,c='orange',ha='center',va='top',
        fontweight='bold')
    plt.draw()
    
# Fill a DataFrame with info to make the TEF extractions
df = pd.DataFrame()
# and plot the stairsteps along the way
s_list = []
i_list = []
j_list = []
ir0_list = []
ir1_list = []
jr0_list = []
jr1_list = []
uv_list = []
pm_list = []
for sn in sn_list:
    plot_line(sn)
    ix, iy = get_stairstep(sn)
    # trim overlap
    ix, iy, trimmed = trim(ix, iy)
    if trimmed > 0:
        print('\n%s: trimmed %d points' % (sn, trimmed))
    # check for duplicates
    xy = list(zip(ix,iy))
    xys = set(xy)
    if len(xy) != len(xys):
        print('\nWarning: there are duplicates in %s' % (sn))
        print(' len(xy) = %d' % (len(xy)))
        print(' len(set(xy)) = %d' % (len(xys)))
    NP = len(ix)
    for ii in range(NP-1):
        ax.plot([lop[ix[ii]],lop[ix[ii+1]]], [lap[iy[ii]],lap[iy[ii+1]]],'-sr',alpha=.5)
        i0 = ix[ii]; i1 = ix[ii+1]
        j0 = iy[ii]; j1 = iy[ii+1]
        if (i0 == i1) and (j1 == j0 + 1): # S to N segment
            if m[j1,i1] + m[j1,i1+1] == 2: # there is ocean on both sides
                i_list.append(i1); j_list.append(j1)
                ir0_list.append(i1); ir1_list.append(i1+1)
                jr0_list.append(j1); jr1_list.append(j1)
                uv_list.append('u')
                pm_list.append(-1)
                s_list.append(sn)
        elif (i0 == i1) and (j1 == j0 - 1): # N to S segment
            if m[j0,i0] + m[j0,i0+1] == 2: # there is ocean on both sides
                i_list.append(i0); j_list.append(j0)
                ir0_list.append(i0); ir1_list.append(i0+1)
                jr0_list.append(j0); jr1_list.append(j0)
                uv_list.append('u')
                pm_list.append(1)
                s_list.append(sn)
        elif (j0 == j1) and (i1 == i0 + 1): # W to E segment
            if m[j1,i1] + m[j1+1,i1] == 2: # there is ocean on both sides
                i_list.append(i1); j_list.append(j1)
                ir0_list.append(i1); ir1_list.append(i1)
                jr0_list.append(j1); jr1_list.append(j1+1)
                uv_list.append('v')
                pm_list.append(1)
                s_list.append(sn)
        elif (j0 == j1) and (i1 == i0 - 1): # E to W segment
            if m[j0,i0] + m[j0+1,i0] == 2: # there is ocean on both sides
                i_list.append(i0); j_list.append(j0)
                ir0_list.append(i0); ir1_list.append(i0)
                jr0_list.append(j0); jr1_list.append(j0+1)
                uv_list.append('v')
                pm_list.append(-1)
                s_list.append(sn)
                
df['sn'] = s_list
df['i'] = i_list
df['j'] = j_list
df['ir0'] = ir0_list
df['jr0'] = jr0_list
df['ir1'] = ir1_list
df['jr1'] = jr1_list
df['uv'] = uv_list
df['pm'] = pm_list

ax.plot(lor[ir0_list],lar[jr0_list],'ob')
ax.plot(lor[ir1_list],lar[jr1_list],'ob')

ax.plot(lou[df.loc[(df.uv=='u') & (df.pm==1),'i'].to_numpy()],
    lau[df.loc[(df.uv=='u') & (df.pm==1),'j'].to_numpy()],'>y')
ax.plot(lou[df.loc[(df.uv=='u') & (df.pm==-1),'i'].to_numpy()],
    lau[df.loc[(df.uv=='u') & (df.pm==-1),'j'].to_numpy()],'<y')
ax.plot(lov[df.loc[(df.uv=='v') & (df.pm==1),'i'].to_numpy()],
    lav[df.loc[(df.uv=='v') & (df.pm==1),'j'].to_numpy()],'^y')
ax.plot(lov[df.loc[(df.uv=='v') & (df.pm==-1),'i'].to_numpy()],
    lav[df.loc[(df.uv=='v') & (df.pm==-1),'j'].to_numpy()],'vy')
    
plt.draw()
    
    
