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

tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'
sect_df = pd.read_pickle(tef2_dir / ('sect_df_' + gctag + '.p'))

if testing:
    sn_list = ['mb8']
    sect_df = sect_df.loc[(sect_df.sn == 'mb7') | 
                          (sect_df.sn == 'mb8') |
                          (sect_df.sn == 'mb9'),:].copy()
    sect_df = sect_df.reset_index(drop=True)
    # I think we get three sections so that we can find the axis limits for plotting.
    # In this case the section we are testing is mb8, and this gives rise to
    # two segments, on on either side, bounded by mb7 to the north and mb9 to the south.
else:
    sn_list = list(sect_df.sn.unique())
    sn_list.remove('jdf1')
    sn_list.remove('sog7')

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
    ax.axis([lor[sect_df.irp.min()-5],lor[sect_df.irp.max()+5],
        lar[sect_df.jrp.min()-5],lar[sect_df.jrp.max()+5]])
    ax.text(.05,.95,gridname,transform=ax.transAxes,
        fontweight='bold')

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
    return ax

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
                print('hit ' + snx + ' on plus side')
                dfx = sect_df.loc[sect_df.sn==snx,:]
                mm[dfx.jrm,dfx.irm] = False
                df_not = df_not.loc[df_not.sn!=snx,:]
                sns_list.append(snx + '_p')
                
            elif len(m) > 0:
                snx = m.sn.values[0]
                print('hit ' + snx + ' on minus side')
                dfx = sect_df.loc[sect_df.sn==snx,:]
                mm[dfx.jrp,dfx.irp] = False
                df_not = df_not.loc[df_not.sn!=snx,:]
                sns_list.append(snx + '_m')
                
        next_ji_list = []
        counter += 1
    print('points = ' + str(len(full_ji_list)))
    return full_ji_list, sns_list

ji_dict = {}    
for sn in sn_list:
    if testing and sn == sn_list[0]:
        ax = initialize_plot()
        
    for pm in [-1, 1]:
        if pm == 1:
            pm_str = 'p'
        elif pm == -1:
            pm_str = 'm'
        sns = sn + '_' + pm_str
        # only do this segment if it is not already done
        done_list = []
        for k in ji_dict.keys():
            done_list += ji_dict[k]['sns_list']
        if sns not in done_list:
            full_ji_list, sns_list = update_mm(sn, pm, m, sect_df)
            ji_dict[sns] = {'ji_list':full_ji_list, 'sns_list': sns_list}
            if testing:
                # plotting to check
                jj = [item[0] for item in full_ji_list]
                ii = [item[1] for item in full_ji_list]
                ax.plot(lor[ii], lar[jj],'sw')
        else:
            print('\nSkipping ' + sns)
if testing:
    plt.show()
#
print(50*'=')
for k in ji_dict.keys():
    print('\n'+k)
    print(ji_dict[k]['sns_list'])
   