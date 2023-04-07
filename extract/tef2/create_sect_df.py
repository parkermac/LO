"""
Create a DataFrame of all the section index and sign info required for section extractions.

Example call:
run create_sect_df -gctag cas6_c0 -test True

"""

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cmocean import cm
from time import time

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# gridname and tag for collection folder
parser.add_argument('-gctag', default='cas6_c0', type=str)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input and output locations
gctag = args.gctag
gridname = gctag.split('_')[0]
Ldir = Lfun.Lstart(gridname=gridname)
grid_fn = Ldir['grid'] / 'grid.nc'
out_name = 'sect_df_' + gctag + '.p'
out_dir = Ldir['LOo'] / 'extract' / 'tef2'
out_fn = out_dir / out_name
collection_dir = out_dir / ('sections_' + gctag)
    
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

def get_sn_list():
    df_list = list(collection_dir.glob('*.p'))
    sn_list = [item.name.replace('.p','') for item in df_list]
    return sn_list

sn_list = get_sn_list()

if args.testing == True:
    sn_list = ['mb8','mb9']
    
# initialize plot
plt.close('all')
if False:
    fig = plt.figure(figsize=(8,8)) # laptop size
else:
    fig = plt.figure(figsize=(12,12)) # external monitor size
ax = fig.add_subplot(111)
ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
pfun.dar(ax)
ax.axis([-127.5,-122,46.5,51])
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
    alpha = .3
    x, y = get_xy(sn)
    ax.plot(x,y,'-o',c='orange', lw=2, alpha=alpha)
    ax.text(x[0],y[0],'\n'+sn,c='orange',ha='center',va='top',
        fontweight='bold')
    plt.draw()
    
# Fill a DataFrame with info to make the TEF extractions
df = pd.DataFrame()

# initialize lists

s_list = [] # section name
i_list = [] # column index of u or v
j_list = [] # row index of u or v

# lists of row and colums indices on the rho grid, with "p" and "m" indicating +/-
irp_list = [] 
irm_list = []
jrp_list = []
jrm_list = []

uv_list = [] # u or v
pm_list = [] # 1 or -1 using right hand rule along section
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
        ax.plot([lop[ix[ii]],lop[ix[ii+1]]], [lap[iy[ii]],lap[iy[ii+1]]],'-sr',alpha=.1)
        i0 = ix[ii]; i1 = ix[ii+1]
        j0 = iy[ii]; j1 = iy[ii+1]
        if (i0 == i1) and (j1 == j0 + 1): # S to N segment
            if m[j1,i1] + m[j1,i1+1] == 2: # there is ocean on both sides
                i_list.append(i1); j_list.append(j1)
                irp_list.append(i1); irm_list.append(i1+1)
                jrp_list.append(j1); jrm_list.append(j1)
                uv_list.append('u')
                pm_list.append(-1)
                s_list.append(sn)
        elif (i0 == i1) and (j1 == j0 - 1): # N to S segment
            if m[j0,i0] + m[j0,i0+1] == 2: # there is ocean on both sides
                i_list.append(i0); j_list.append(j0)
                irm_list.append(i0); irp_list.append(i0+1)
                jrm_list.append(j0); jrp_list.append(j0)
                uv_list.append('u')
                pm_list.append(1)
                s_list.append(sn)
        elif (j0 == j1) and (i1 == i0 + 1): # W to E segment
            if m[j1,i1] + m[j1+1,i1] == 2: # there is ocean on both sides
                i_list.append(i1); j_list.append(j1)
                irm_list.append(i1); irp_list.append(i1)
                jrm_list.append(j1); jrp_list.append(j1+1)
                uv_list.append('v')
                pm_list.append(1)
                s_list.append(sn)
        elif (j0 == j1) and (i1 == i0 - 1): # E to W segment
            if m[j0,i0] + m[j0+1,i0] == 2: # there is ocean on both sides
                i_list.append(i0); j_list.append(j0)
                irp_list.append(i0); irm_list.append(i0)
                jrp_list.append(j0); jrm_list.append(j0+1)
                uv_list.append('v')
                pm_list.append(-1)
                s_list.append(sn)
                
df['sn'] = s_list
df['i'] = i_list
df['j'] = j_list
df['irp'] = irp_list
df['jrp'] = jrp_list
df['irm'] = irm_list
df['jrm'] = jrm_list
df['uv'] = uv_list
df['pm'] = pm_list

# save pickled DataFrame
df.to_pickle(out_fn)

# plot to check results
dfc = pd.read_pickle(out_fn)

ax.plot(lor[dfc.irp],lar[dfc.jrp],'or')
ax.plot(lor[dfc.irm],lar[dfc.jrm],'ob')

ax.plot(lou[dfc.loc[(dfc.uv=='u') & (dfc.pm==1),'i']],
    lau[dfc.loc[(dfc.uv=='u') & (dfc.pm==1),'j']],'>y')
ax.plot(lou[dfc.loc[(dfc.uv=='u') & (dfc.pm==-1),'i']],
    lau[dfc.loc[(dfc.uv=='u') & (dfc.pm==-1),'j']],'<y')
ax.plot(lov[dfc.loc[(dfc.uv=='v') & (dfc.pm==1),'i']],
    lav[dfc.loc[(dfc.uv=='v') & (dfc.pm==1),'j']],'^y')
ax.plot(lov[dfc.loc[(dfc.uv=='v') & (dfc.pm==-1),'i']],
    lav[dfc.loc[(dfc.uv=='v') & (dfc.pm==-1),'j']],'vy')
    
plt.draw()
    
    
