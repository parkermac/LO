"""
Plot results of a particle tracking experiment, specific to experiments about
the residence time of Hood Canal.  Integrates particle tracking results with
the cas6_v3_lo8dye run.
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import pickle
import pandas as pd
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.path as mpltPath
from time import time
from datetime import datetime, timedelta

Ldir = Lfun.Lstart()

# get release Dataset
indir0 = Ldir['LOo'] + 'tracks/'
indir = 'AllHC3d_ndiv12_3d_3mo/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)

NT, NP = dsr['lon'].shape

# get a list of datetimes
ot_vec = dsr['ot'][:]
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]
t = (ot_vec - ot_vec[0])/3600

# Gather particle data
# packed [time, particle #]
x = dsr['lon'][:]
y = dsr['lat'][:]
#z = dsr['z'][:]
h = dsr['h'][:]
u = dsr['u'][:]
v = dsr['v'][:]
w = dsr['w'][:]
salt = dsr['salt'][:]
temp = dsr['temp'][:]
cs = dsr['cs'][:]
zeta = dsr['zeta'][:]
dsr.close()

# this version of z ignores the tide and always has z <= 0
z = cs*h

dsg = nc4.Dataset(indir0 + indir + 'grid.nc')
lon = dsg['lon_rho'][:].data
lat = dsg['lat_rho'][:].data
mask = dsg['mask_rho'][:].data # 1 = water
dsg.close()

# load the TEF segment ji_dicts and volumes
voldir = Ldir['LOo'] + 'tef2/volumes_cas6/'
v_df = pd.read_pickle(voldir + 'volumes.p')
bathy_dict = pickle.load(open(voldir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(voldir + 'ji_dict.p', 'rb'))
seg_list = list(v_df.index)
hc_segs = ['H' + str(item) for item in range(1,9)]

def packer(x,y):
    # creates 2D array out of vectors, with columns x,y
    N = len(x)
    xy = np.concatenate((x.reshape(N,1), y.reshape(N,1)), axis=1)
    return xy

# Then, count all the particles in all segments for each hour
# and save the results in a DataFrame OR load that DataFrame.
count_fn = indir0 + indir + 'np_df.p'
if os.path.isfile(count_fn):
    np_df = pd.read_pickle(count_fn)
else:
    # This takes 30-160 sec per month on my mac
    # ==============================================================

    j_dict = {}; i_dict = {}
    for seg_name in seg_list:
        jj = []; ii = []
        ji_list_full = ji_dict[seg_name]
        for ji in ji_list_full:
            jj.append(ji[0])
            ii.append(ji[1])
        jjj = np.array(jj)
        iii = np.array(ii)
        j_dict[seg_name] = jjj
        i_dict[seg_name] = iii
    
    # now find the convex hull around the segments, and the points in that polygon
    tt0 = time()
    path_dict = {}
    for sn in seg_list:
        xx = lon[j_dict[sn], i_dict[sn]]
        yy = lat[j_dict[sn], i_dict[sn]]
        nseg = len(xx)
        points = packer(xx,yy)
        hull = ConvexHull(points)
        # these xy, yh are the points that define the polygon of the hull
        xh = points[hull.vertices,0]
        yh = points[hull.vertices,1]
        xyh = packer(xh, yh)
        path_dict[sn] = mpltPath.Path(xyh)
    print('Time to get all hulls = %0.2f sec' % (time()-tt0))
    
    # and count particles in each segment at each time (longest step)
    np_df = pd.DataFrame(index=dt_list, columns = seg_list)
    tt0 = time()
    for tt in range(NT):
        if np.remainder(tt,24) == 0:
            print(' - working on %s' % (str(dt_list[tt])))
        x1 = x[tt,:]
        y1 = y[tt,:]
        xyp = packer(x1, y1)
        for sn in seg_list:
            path = path_dict[sn]
            is_inside = path.contains_points(xyp)
            # x11 = x1[is_inside]
            # y11 = y1[is_inside]
            np_df.loc[dt_list[tt],sn] = is_inside.sum()
            #print('  - %s has %d points (%d points)' % (sn, len(x11), is_inside.sum()))
    print('Time to count points in all segments = %0.2f sec' % (time()-tt0))
    np_df.to_pickle(count_fn)

# ==============================================================

# manipulate particle tallies
np_hc_df = np_df[hc_segs].copy()
np_norm_hc_df = np_hc_df/(np_hc_df.sum(axis=1)[0])
np_HC_df = np_hc_df.sum(axis=1)
np_norm_HC_df = np_HC_df/np_HC_df[0]

# also load and manipulate the dye results
dye_indir = Ldir['LOo'] + 'tef2/cas6_v3_lo8dye_2019.07.04_2019.10.04/flux/'
dye_df = pd.read_pickle(dye_indir + 'hourly_segment_dye.p') # mean concentration in segments
dye_vol_df = pd.read_pickle(dye_indir + 'hourly_segment_dye_volume.p')
dye_hc_df = dye_df[hc_segs].copy()
dye_vol_hc_df = dye_vol_df[hc_segs].copy()
#
net_dye_hc_df = dye_hc_df * dye_vol_hc_df
dye_HC_df = net_dye_hc_df.sum(axis=1)/dye_vol_hc_df.sum(axis=1)
#
dye_norm_hc_df = net_dye_hc_df/net_dye_hc_df.sum(axis=1)[0]

# and load a comparable flux_engine result
cc_indir = Ldir['LOo'] + 'tef/flux_engine/cas6_v3_lo8b/'
cc = pd.read_pickle(cc_indir + 'IC_HoodCanal_2019_Summer.p')
#cc = pd.read_pickle(cc_indir + 'IC_HoodCanal_2019_Fall.p')
cc = cc.loc[:92,:]
cil = cc.index.to_list()
cct = [(datetime(2019,7,4)+timedelta(days=item)) for item in cil]
CC = pd.DataFrame(index=cct, columns=hc_segs)
cc_new = pd.DataFrame(cc.values, index=cct, columns=cc.columns)
v_ser = v_df.loc[hc_segs,'volume m3']
V = v_ser.sum()
for sn in hc_segs:
    v = v_ser[sn]
    frac = .2
    v_f = frac*v
    v_s = (1-frac)*v
    CC.loc[:,sn] = v_f*cc_new.loc[:,sn+'_f'] + v_s*cc_new.loc[:,sn+'_s']
CC = CC/V
CC_net = CC.sum(axis=1)

# PLOTTING
plt.close('all')
fs = 16
plt.rc('font', size=fs)
fig = plt.figure(figsize=(18,10))

# Map
ax = fig.add_subplot(121)
ax.plot(x[NT-1,:],y[NT-1,:], '.', color='k', ms=1, alpha=.1)
# the lines below were for testing the path.contains_points(xyp) call
# RESULT: it works great!
#ax.plot(xh, yh, 'r-', lw=2)
#ax.plot(x11,y11, '.', color='r', ms=2, alpha=.6)
pfun.dar(ax)
pfun.add_coast(ax)
aa = [-124, -122, 47, 49]
ax.axis(aa)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks([-124, -123, -122])
ax.set_yticks([47, 48, 49])
ax.set_title(indir.replace('/',''))

c_list = ['r','b','g','y','m','c','k','orange']
c_dict = dict(zip(hc_segs,c_list))
for sn in hc_segs:
    ax.text(v_df.loc[sn,'lon'],v_df.loc[sn,'lat'],sn,c=c_dict[sn],weight='bold', ha='center', va='center')
ax.text(.95,.9,'(a)',weight='bold',transform=ax.transAxes, ha='right')

# time series
ax = fig.add_subplot(222)
lw=1
np_norm_HC_df.plot(ax=ax, lw=lw, color='k', ls='--', label='Particles', legend=True)
dye_HC_df.plot(ax=ax, lw=lw, color='k', ls='-', label='Dye', legend=True)
CC_net.plot(ax=ax, lw=lw, color='k', ls=':', label='Box Model', legend=True)
ax.set_ylim(0,1)
ax.set_xticklabels([])
ax.text(.95,.9,'(b)',weight='bold',transform=ax.transAxes, ha='right')

ax = fig.add_subplot(224)
np_norm_hc_df.plot(ax=ax, ls='--', color=c_list, alpha=1, legend=False)
dye_norm_hc_df.plot(ax=ax, ls='-', color=c_list, alpha=1, legend=False)
CC.plot(ax=ax, ls=':', color=c_list, alpha=1, legend=False)
ax.text(.95,.9,'(c)',weight='bold',transform=ax.transAxes, ha='right')


plt.show()
plt.rcdefaults()


