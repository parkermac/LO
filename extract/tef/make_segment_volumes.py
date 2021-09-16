"""
Find the volume of each segment
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import netCDF4 as nc

from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun
import tef_fun
import flux_fun

gridname = 'cas6'
Ldir = Lfun.Lstart(gridname=gridname)

# output location
outdir = Ldir['LOo'] / 'extract' / 'tef' / ('volumes_' + Ldir['gridname'])
Lfun.make_dir(outdir, clean=True)

fng = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(fng, only_G=True)
h = np.ma.masked_where(G['mask_rho']==0, G['h'])
x = G['lon_rho']
y = G['lat_rho']
xp, yp = pfun.get_plon_plat(x,y)
m = G['mask_rho'] == 1 # Boolean, True = water
DA = G['DX'] * G['DY']

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df(gridname)

testing = False

# segment definitions, assembled by looking at the figure
# created by flux_seg_map.py
segs = flux_fun.segs

if testing == True:
    seg_name_list = ['J4']
    plt.close('all')
    # start a useful plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.pcolormesh(xp, yp, h, cmap='terrain_r', vmin=-100, vmax = 400, alpha=.1)
    pfun.dar(ax)
    pfun.add_coast(ax)
else:
     seg_name_list = segs.keys()

# initialize a DataFrame to hold all volumes:
vol_df = pd.DataFrame(index=seg_name_list, columns=['volume m3', 'area m2', 'lon', 'lat'])

# create a dict to hold the bathy file
bathy_dict = {'h':h, 'xp':xp, 'yp':yp}

# initialize a dict to store all the ji lists that define each segment area
ji_dict = {}

for seg_name in seg_name_list:

    print('Segment: ' + seg_name)

    # get the names of the sections around this segment
    seg = segs[seg_name]

    # initialize a DataFrame to hold segment info
    seg_df = pd.DataFrame(columns=['ii0','ii1','jj0','jj1',
        'sdir','side','Lon0','Lon1','Lat0','Lat1',])

    # fill the DataFrame for this segment
    for side in list('SNWE'):
        sect_list = seg[side]
        for sect_name in sect_list:
            # get section lat, lon, and other info
            x0, x1, y0, y1 = sect_df.loc[sect_name,:]
            # get indices for this section
            ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
            Lon0 = Lon.min(); Lon1 = Lon.max()
            Lat0 = Lat.min(); Lat1 = Lat.max()
            seg_df.loc[sect_name,'ii0'] = ii0
            seg_df.loc[sect_name,'ii1'] = ii1
            seg_df.loc[sect_name,'jj0'] = jj0
            seg_df.loc[sect_name,'jj1'] = jj1
            seg_df.loc[sect_name,'sdir'] = sdir
            seg_df.loc[sect_name,'side'] = side
            seg_df.loc[sect_name,'Lon0'] = Lon0
            seg_df.loc[sect_name,'Lon1'] = Lon1
            seg_df.loc[sect_name,'Lat0'] = Lat0
            seg_df.loc[sect_name,'Lat1'] = Lat1


    if testing:
        # focus the plot axes around this segment
        pad = .5
        ax.axis([seg_df['Lon0'].min()-pad,seg_df['Lon1'].max()+pad,
            seg_df['Lat0'].min()-pad,seg_df['Lat1'].max()+pad])

    # initialize a mask
    mm = m.copy() # boolean array, True over water

    # initialize some lists
    full_ji_list = [] # full list of indices of good rho points inside the volume
    this_ji_list = [] # current list of indices of good rho points inside the volume
    next_ji_list = [] # next list of indices of good rho points inside the volume

    for sn in seg_df.index:
        s = seg_df.loc[sn,:]
        # mask the rho grid points on the outside of the TEF sections
        if s['sdir'] == 'NS' and s['side'] == 'W':
            mm[s['jj0']:s['jj1']+1, s['ii0']] = False
        elif s['sdir'] == 'NS' and s['side'] == 'E':
            mm[s['jj0']:s['jj1']+1, s['ii1']] = False
        if s['sdir'] == 'EW' and s['side'] == 'S':
            mm[s['jj0'], s['ii0']:s['ii1']+1] = False
        if s['sdir'] == 'EW' and s['side'] == 'N':
            mm[s['jj1'], s['ii0']:s['ii1']+1] = False
        # doing this will form a natural barrier for the "search robot"
#
    if testing:
        # same as the loop above, but for plottting
        for sn in seg_df.index:
            s = seg_df.loc[sn,:]
            # the black dots show rho gridpoints OUTSIDE of the segment volume, while
            # the magenta ones are those just inside
            if s['sdir'] == 'NS' and s['side'] == 'W':
                ax.plot(x[s['jj0']:s['jj1']+1, s['ii0']], y[s['jj0']:s['jj1']+1, s['ii0']], 'ok')
                ax.plot(x[s['jj0']:s['jj1']+1, s['ii1']], y[s['jj0']:s['jj1']+1, s['ii1']], 'om')
            elif s['sdir'] == 'NS' and s['side'] == 'E':
                ax.plot(x[s['jj0']:s['jj1']+1, s['ii0']], y[s['jj0']:s['jj1']+1, s['ii0']], 'om')
                ax.plot(x[s['jj0']:s['jj1']+1, s['ii1']], y[s['jj0']:s['jj1']+1, s['ii1']], 'ok')
            if s['sdir'] == 'EW' and s['side'] == 'S':
                ax.plot(x[s['jj0'], s['ii0']:s['ii1']+1], y[s['jj0'], s['ii0']:s['ii1']+1], 'ok')
                ax.plot(x[s['jj1'], s['ii0']:s['ii1']+1], y[s['jj1'], s['ii0']:s['ii1']+1], 'om')
            if s['sdir'] == 'EW' and s['side'] == 'N':
                ax.plot(x[s['jj0'], s['ii0']:s['ii1']+1], y[s['jj0'], s['ii0']:s['ii1']+1], 'om')
                ax.plot(x[s['jj1'], s['ii0']:s['ii1']+1], y[s['jj1'], s['ii0']:s['ii1']+1], 'ok')
#
    # deploy the "search robot" flux_fun.update_mm()
    # the algorithm is that we start the search at a good rho point at either end of each TEF
    # section, and allow it to fill in as much of the remaining points as it can
    for sn in seg_df.index:
        #print(sn)
        s = seg_df.loc[sn,:]
        if s['sdir'] == 'NS' and s['side'] == 'W':
            ji = (s['jj0'],s['ii1'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
            ji = (s['jj1'],s['ii1'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
        elif s['sdir'] == 'NS' and s['side'] == 'E':
            ji = (s['jj0'],s['ii0'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
            ji = (s['jj1'],s['ii0'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
        elif s['sdir'] == 'EW' and s['side'] == 'S':
            ji = (s['jj1'],s['ii0'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
            ji = (s['jj1'],s['ii1'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
        elif s['sdir'] == 'EW' and s['side'] == 'N':
            ji = (s['jj0'],s['ii0'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
            ji = (s['jj0'],s['ii1'])
            mm, this_ji_list, full_ji_list, next_ji_list = flux_fun.update_mm(ji, mm,
                    this_ji_list, full_ji_list, next_ji_list)
    # check results for extras
    if len(set(full_ji_list)) != len(full_ji_list):
        print(' -- Warning: had to remove duplicates from list')
        full_ji_list = set(full_ji_list)
        
    ji_dict[seg_name] = full_ji_list

    if testing:
        hh = np.zeros(h.shape)
        for ji in full_ji_list: 
            hh[ji]=1
        H = hh==0
        hm = np.ma.masked_where(H,h)
        ax.pcolormesh(xp, yp, hm, cmap='terrain_r', vmin=-100, vmax = 400)
        
#
    # find the volume and surface area
    volume = 0
    area = 0
    lon = 0
    lat = 0
    for ji in full_ji_list:
        area += DA[ji]
        volume += h[ji] * DA[ji]
        lon += x[ji]
        lat += y[ji]
    lon = lon / len(full_ji_list)
    lat = lat / len(full_ji_list)
    print(' -- Area = %0.1f km2' % (area/1e6))
    print(' -- Volume = %0.1f km3' % (volume/1e9))

    vol_df.loc[seg_name,'volume m3'] = volume
    vol_df.loc[seg_name,'area m2'] = area
    vol_df.loc[seg_name,'lon'] = lon
    vol_df.loc[seg_name,'lat'] = lat
    
# repack the ji info into separate dicts: keys are segment names, and values are
# 1-D arrays of the indices, ready to use for fancy indexing.
j_dict = {}; i_dict = {}
for seg_name in seg_name_list:
    jj = []; ii = []
    ji_list_full = ji_dict[seg_name]
    for ji in ji_list_full:
        jj.append(ji[0])
        ii.append(ji[1])
    jjj = np.array(jj)
    iii = np.array(ii)
    j_dict[seg_name] = jjj
    i_dict[seg_name] = iii
    
if not testing:
    vol_df.to_pickle(outdir / 'volumes.p')
    pickle.dump(bathy_dict, open(outdir / 'bathy_dict.p', 'wb'))
    pickle.dump(i_dict, open(outdir / 'i_dict.p', 'wb'))
    pickle.dump(j_dict, open(outdir / 'j_dict.p', 'wb'))
#
if testing:
    plt.show()


