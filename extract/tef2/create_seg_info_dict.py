"""
Code to define tef2 segments. Specifically we get all the j,i indices on
the rho grid for all the regions between sections. for each segment we also
get a list of the bounding sections and the rivers.

The results are saved in a dict of dicts called "seg_info_dict":

seg_info_dict has one key for each segment, e.g. 'mb8_m'
then seg_info_dict['mb8_m'] is a dict with three keys:
['ji_list', 'sns_list', 'riv_list']
 - ji_list is a list of tuples (j,i) that are indices of points in the rho grid in the segment
 - sns_list is a list of all the bounding sections, including their signs. Note that the segment
   key is one member of this list
 - riv_list is a list of the rivers or point sources entering the segment

It also creates and saves a DataFrame with columns:
['volume m3', 'area m2', 'lon', 'lat']
which are the volume, surface area, and mean lon and lat of each segment.

To test on mac:

run create_seg_info_dict -gctag cas6_c0 -riv riv00 -small True -test True

or just:

run create_seg_info_dict (uses all defaults)

Performance: takes a few minutes to run for real on cas6_c0_traps2. Or just 30 sec
for _riv00.

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
import pickle

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# gridname and tag for collection folder
parser.add_argument('-gctag', default='cas6_c0', type=str)
parser.add_argument('-riv', default='riv00', type=str)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input and output locations
gctag = args.gctag
gridname = gctag.split('_')[0]
riv = args.riv
testing = args.testing
Ldir = Lfun.Lstart(gridname=gridname)
grid_fn = Ldir['grid'] / 'grid.nc'
if len(riv) > 0:
    out_name = 'seg_info_dict_' + gctag + '_' + riv + '.p'
else:
    out_name = 'seg_info_dict_' + gctag + '_noriv.p'
out_dir = Ldir['LOo'] / 'extract' / 'tef2'
Lfun.make_dir(out_dir)
out_fn = out_dir / out_name
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
# these are used for making segment volumes
H = h.copy()
DX = 1/ds.pm.values
DY = 1/ds.pn.values
DA = DX * DY
DV = H * DA
lon_rho = ds.lon_rho.values
lat_rho = ds.lat_rho.values
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

sect_df = pd.read_pickle(out_dir / ('sect_df_' + gctag + '.p'))

# Get river info if there is any. This must have been created from a rivers.nc
# file using create_river_info.py.
if len(riv) > 0:
    riv_name = 'riv_df_' + gridname + '_' + riv + '.p'
    riv_fn = out_dir / riv_name
    riv_df = pd.read_pickle(riv_fn)
    riv_names = list(riv_df.name)
    riv_i = riv_df.irho.to_numpy(dtype = int)
    riv_j = riv_df.jrho.to_numpy(dtype = int)
    riv_ji = [(riv_j[ii], riv_i[ii]) for ii in range(len(riv_i))]
    # make two dicts out of this info, just reversing the key-value order
    riv_ji_dict = dict(zip(riv_ji, riv_names)) # used to search for river names from ji points
    riv_ji_dict2 = dict(zip(riv_names, riv_ji)) # used to plot rivers from names
    do_riv = True
else:
    print('Not using river info.')
    do_riv = False

if testing:
    if False:
        sn_list = ['mb8','mb9']
        sect_df = sect_df.loc[(sect_df.sn == 'mb7') | 
                              (sect_df.sn == 'mb8') |
                              (sect_df.sn == 'mb9') |
                              (sect_df.sn == 'mb10') |
                              (sect_df.sn == 'tn1'),:].copy()
    else:
        sn_list = ['jdf2']
        sect_df = sect_df.loc[(sect_df.sn == 'jdf1') | 
                              (sect_df.sn == 'jdf2') |
                              (sect_df.sn == 'jdf3'),:].copy()
        
    sect_df = sect_df.reset_index(drop=True)
    # We need to have all neighboring sections in sect_df for plotting.
else:
    sn_list = list(sect_df.sn.unique())

# NOTE: We have to remove the open boundary sections from the list
# or else our algorithm will make a segment out of the ocean!
bounding_section_fn = out_dir / ('sections_' + gctag) / 'bounding_sections.txt'
with open(bounding_section_fn,'r') as f:
    bounding_section_list = f.read().split('\n')
for sn in bounding_section_list:
    try:
        sn_list.remove(sn)
    except ValueError:
        pass

# start of routine to find all rho-grid points within a segment

def update_mm(sn, pm, m, sect_df):
    
    # initialize some lists
    full_ji_list = [] # full list of indices of good rho points inside the volume
    this_ji_list = [] # current list of indices of good rho points inside the volume
    next_ji_list = [] # next list of indices of good rho points inside the volume
    sns_list = [] # list of bounding sections and signs of interior
    
    if pm == 1:
        pm_str = 'p'
    elif pm == -1:
        pm_str = 'm'
    print('\nSection %s, %s side' % (sn, pm_str))
    
    sns_list.append(sn + '_' + pm_str)
    
    df = sect_df.loc[sect_df.sn==sn,:]
    
    df_not = sect_df.loc[sect_df.sn!=sn,:]
    
    # initialize a mask
    mm = m == 1 # boolean array, True over water
    if pm == 1:
        mm[df.jrm,df.irm] = False
        ji = (df.loc[df.index[0],'jrp'],df.loc[df.index[0],'irp'])
    elif pm == -1:
        mm[df.jrp,df.irp] = False
        ji = (df.loc[df.index[0],'jrm'],df.loc[df.index[0],'irm'])
    
    if mm[ji] == True:
        this_ji_list.append(ji)
        full_ji_list.append(ji)
    for ji in this_ji_list:
        mm[ji] = False
    
    counter = 0
    while len(this_ji_list) > 0:
        # print('iteration ' + str(counter))
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
                print(' - hit ' + snx + ' on plus side')
                dfx = sect_df.loc[sect_df.sn==snx,:]
                mm[dfx.jrm,dfx.irm] = False
                df_not = df_not.loc[df_not.sn!=snx,:]
                sns_list.append(snx + '_p')
                
            elif len(m) > 0:
                snx = m.sn.values[0]
                print(' - hit ' + snx + ' on minus side')
                dfx = sect_df.loc[sect_df.sn==snx,:]
                mm[dfx.jrp,dfx.irp] = False
                df_not = df_not.loc[df_not.sn!=snx,:]
                sns_list.append(snx + '_m')
                
        next_ji_list = []
        counter += 1
    print(' - points = ' + str(len(full_ji_list)))
    return full_ji_list, sns_list
    
def find_riv_list(ji_list):
    riv_list = []
    for ji in ji_list:
        if ji in riv_ji_dict.keys():
            riv_list.append(riv_ji_dict[ji])
    return riv_list

# this it the loop where we actually find the segment info
tt0 = time()
seg_info_dict = {}
for sn in sn_list:
        
    for pm in [-1, 1]:
        if pm == 1:
            pm_str = 'p'
        elif pm == -1:
            pm_str = 'm'
        sns = sn + '_' + pm_str
        
        # only do this segment if it is not already done
        done_list = [] # we regenerate this list for each attempt
        for k in seg_info_dict.keys():
            done_list += seg_info_dict[k]['sns_list']
        if sns not in done_list:
            full_ji_list, sns_list = update_mm(sn, pm, m, sect_df)
            seg_info_dict[sns] = {'ji_list':full_ji_list, 'sns_list': sns_list}
            
            # find rivers in this segment
            riv_list = []
            if do_riv == True:
                riv_list = find_riv_list(full_ji_list)
            seg_info_dict[sns]['riv_list'] = riv_list
            
        else:
            print('\nSkipping ' + sns)
if not testing:
    pickle.dump(seg_info_dict, open(out_fn,'wb'))
print('\nElapsed time = %0.1f seconds' % (time()-tt0))
print('output file: %s' % (str(out_fn)))

# calculate volumes, areas, and mean lon,lat of all segments
# initialize a DataFrame to hold all volumes:
sns_list = seg_info_dict.keys()
vol_df = pd.DataFrame(index=sns_list, columns=['volume m3', 'area m2', 'lon', 'lat'])
for sns in sns_list:
    ji_list = seg_info_dict[sns]['ji_list']
    # make index vectors for fancy indexing
    jj = []; ii = []
    for ji in ji_list:
        jj.append(ji[0])
        ii.append(ji[1])
    JJ = np.array(jj,dtype=int)
    II = np.array(ii,dtype=int)
    volume = DV[JJ,II].sum()
    area = DA[JJ,II].sum()
    lon = lon_rho[JJ,II].mean()
    lat = lat_rho[JJ,II].mean()
    vol_df.loc[sns,'volume m3'] = volume
    vol_df.loc[sns,'area m2'] = area
    vol_df.loc[sns,'lon'] = lon
    vol_df.loc[sns,'lat'] = lat
vol_df.to_pickle(out_dir / ('vol_df_' + gctag + '.p'))

if testing:
    print('\n' + 50*'=')
    for k in seg_info_dict.keys():
        print('\nSegment key: ' + str(k))
        print(' - Bounding sections: ' + str(seg_info_dict[k]['sns_list']))
        if do_riv:
            print(' - Included rivers:')
            this_riv_list = seg_info_dict[k]['riv_list']
            for rn in this_riv_list:
                print('    ' + rn)

def initialize_plot():
    # initialize plot
    plt.close('all')
    if args.small:
        fig = plt.figure(figsize=(8,8)) # laptop size
    else:
        fig = plt.figure(figsize=(12,12)) # external monitor size
    ax = fig.add_subplot(111)
    ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
    pfun.dar(ax)
            
    ax.text(.05,.9,gctag + '\n' + riv, transform=ax.transAxes,
        fontweight='bold')

    ms = 4
    ax.plot(lor[sect_df.irp],lar[sect_df.jrp],'oc',markersize=ms)
    ax.plot(lor[sect_df.irm],lar[sect_df.jrm],'or',markersize=ms)

    ax.plot(lou[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==1),'i']],
        lau[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==1),'j']],'>y',markersize=ms)
    ax.plot(lou[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==-1),'i']],
        lau[sect_df.loc[(sect_df.uv=='u') & (sect_df.pm==-1),'j']],'<y',markersize=ms)
    ax.plot(lov[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==1),'i']],
        lav[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==1),'j']],'^y')
    ax.plot(lov[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==-1),'i']],
        lav[sect_df.loc[(sect_df.uv=='v') & (sect_df.pm==-1),'j']],'vy',markersize=ms)
    
    for sn in sect_df.sn.unique():
        df = sect_df.loc[sect_df.sn==sn,['i','j']].copy()
        df = df.reset_index(drop=True)
        i = df.loc[0,'i']
        j = df.loc[0,'j']
        ax.text(lor[i],lar[j],'\n'+sn,c='orange',ha='center',va='top',
            fontweight='bold',bbox=pfun.bbox)
    return ax
    
if testing:
    ax = initialize_plot()

    # plotting to check
    jjj = []
    iii = []
    for sns in seg_info_dict.keys():
        ji_list = seg_info_dict[sns]['ji_list']
        jj = [item[0] for item in ji_list]
        ii = [item[1] for item in ji_list]
        jjj += jj
        iii += ii
        ax.plot(lor[ii], lar[jj],'sw',markersize=1)


        # add rivers
        if do_riv:
            riv_list = seg_info_dict[sns]['riv_list']
            for rn in riv_list:
                rj, ri = riv_ji_dict2[rn]
                ax.plot(lor[ri], lar[rj], 'om')
               
    imin = np.min(np.array(iii))
    imax = np.max(np.array(iii))
    jmin = np.min(np.array(jjj))
    jmax = np.max(np.array(jjj))
    ax.axis([lor[imin-5],lor[imax+5],
        lar[jmin-5],lar[jmax+5]])
    
   
    plt.show()
   