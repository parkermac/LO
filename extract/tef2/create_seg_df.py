"""
Code to define tef2 segments. Specifically we get all the j,i indices on
the rho grid for all the regions between sections.

To test on mac:
run create_seg_df -gctag cas6_c0 -small True -test True
run create_seg_df -gctag cas6_c0 -test True

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
# gridname and tag for collection folder
parser.add_argument('-gctag', default='cas6_c0', type=str)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input and output locations
gctag = args.gctag
gridname = gctag.split('_')[0]
testing = args.testing
Ldir = Lfun.Lstart(gridname=gridname)
grid_fn = Ldir['grid'] / 'grid.nc'
out_name = 'seg_df_' + gctag + '.p'
out_dir = Ldir['LOo'] / 'extract' / 'tef2'
out_fn = out_dir / out_name
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
h[m==0] = np.nan
# coordinates for plotting
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
# coordinates for convenience in plotting
lop = ds.lon_psi[0,:].values
lap = ds.lat_psi[:,0].values
lor = ds.lon_rho[0,:].values
lar = ds.lat_rho[:,0].values
# lon_rho = ds.lon_rho.values
# lat_rho = ds.lat_rho.values
lou = ds.lon_u[0,:].values
lau = ds.lat_u[:,0].values
lov = ds.lon_v[0,:].values
lav = ds.lat_v[:,0].values
ds.close

tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'
sect_df = pd.read_pickle(tef2_dir / ('sect_df_' + gctag + '.p'))

if testing:
    sect_df = sect_df.loc[(sect_df.sn == 'mb8') | (sect_df.sn == 'mb9'),:].copy()
    sect_df = sect_df.reset_index(drop=True)
    
# initialize plot
plt.close('all')
if args.small:
    fig = plt.figure(figsize=(8,8)) # laptop size
else:
    fig = plt.figure(figsize=(12,12)) # external monitor size
ax = fig.add_subplot(111)
ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
pfun.dar(ax)
ax.axis([lor[sect_df.irp.min()-5],lor[sect_df.irp.max()+5],
    lar[sect_df.jrp.min()-5],lar[sect_df.jrp.max()+5]])
ax.text(.05,.95,gridname,transform=ax.transAxes,
    fontweight='bold')
plt.show()

ax.plot(lor[sect_df.irp],lar[sect_df.jrp],'or')
ax.plot(lor[sect_df.irm],lar[sect_df.jrm],'ob')

ax.plot(lou[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==1),'i']],
    lau[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==1),'j']],'>y')
ax.plot(lou[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==-1),'i']],
    lau[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==-1),'j']],'<y')
ax.plot(lov[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==1),'i']],
    lav[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==1),'j']],'^y')
ax.plot(lov[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==-1),'i']],
    lav[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==-1),'j']],'vy')
    
for sn in sect_df.sn.unique():
    df = sect_df.loc[sect_df.sn==sn,['i','j']].copy()
    df = df.reset_index(drop=True)
    i = df.loc[0,'i']
    j = df.loc[0,'j']
    ax.text(lor[i],lar[j],'\n'+sn,c='orange',ha='center',va='top',
        fontweight='bold')
    
plt.draw()

# start of routine to find all rho-grid within a segment

def update_mm(ji, mm, this_ji_list, full_ji_list, next_ji_list, df_not):
    if mm[ji] == True:
        this_ji_list.append(ji)
        full_ji_list.append(ji)
    for ji in this_ji_list:
        mm[ji] = False
    counter = 0
    while len(this_ji_list) > 0:
        print('iteration ' + str(counter))
        for ji in this_ji_list:
            JI = (ji[0]+1, ji[1]) # North
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0], ji[1]+1) # East
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0]-1, ji[1]) # South
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0], ji[1]-1) # West
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
        for ji in next_ji_list:
            full_ji_list.append(ji)
        this_ji_list = next_ji_list.copy()
        
        # check to see if we have hit another section
        for ji in this_ji_list:
            p = df_not.loc[(df_not.irp==ji[1]) & (df_not.jrp==ji[0]),:]
            m = df_not.loc[(df_not.irm==ji[1]) & (df_not.jrm==ji[0]),:]
            if len(p) > 0:
                snx = p.sn.values[0]
                print('hit ' + snx + ' on plus side')
                dfx = sect_df.loc[sect_df.sn==snx,:].copy()
                mm[dfx.jrm,dfx.irm] = False
                df_not = df_not.loc[df_not.sn!=snx,:].copy()
            elif len(m) > 0:
                snx = m.sn.values[0]
                print('hit ' + snx + ' on minus side')
                dfx = sect_df.loc[sect_df.sn==snx,:].copy()
                mm[dfx.jrp,dfx.irp] = False
                df_not = df_not.loc[df_not.sn!=snx,:].copy()
        
        next_ji_list = []
        counter += 1
    return mm, this_ji_list, full_ji_list, next_ji_list

# initialize a mask
mm = m == 1 # boolean array, True over water

# initialize some lists
full_ji_list = [] # full list of indices of good rho points inside the volume
this_ji_list = [] # current list of indices of good rho points inside the volume
next_ji_list = [] # next list of indices of good rho points inside the volume

if testing:
    sn = 'mb9'
    pm = 1
    df = sect_df.loc[sect_df.sn==sn,:].copy()
    df_not = sect_df.loc[sect_df.sn!=sn,:].copy()
    if pm == 1:
        mm[df.jrm,df.irm] = False
        ji = (df.loc[df.index[0],'jrp'],df.loc[df.index[0],'irp'])
        mm0, this_ji_list0, full_ji_list0, next_ji_list0 = \
            update_mm(ji, mm, this_ji_list, full_ji_list, next_ji_list, df_not)
            
jj = [item[0] for item in full_ji_list]
ii = [item[1] for item in full_ji_list]

ax.plot(lor[ii], lar[jj],'sw')
   