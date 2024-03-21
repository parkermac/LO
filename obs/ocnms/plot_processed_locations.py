# this is a script to plot processed OCNMS data from LO/obs 
# This script will generate a data availabilty plot for 10 sites 

# imports
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

#from time import time
#import sys 
import xarray as xr
import numpy as np
#import netCDF4 as nc

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#import matplotlib.dates as mdates
#from datetime import datetime
import datetime as dt
import pandas as pd

Ldir = Lfun.Lstart()

# processed data location
source = 'ocnms'
otype = 'moor' 
data_dir = Ldir['LOo'] / 'obs' / source / otype

# named in this order for plotting 
sn_name_dict = {
    'MB042':'Makah Bay 42m',
    'MB015':'Makah Bay 15m',
    'CA042':'Cape Alava 42m',
    'CA015':'Cape Alava 15m',
    'TH042':'Teahwhit Head 42m',
    'TH015':'Teahwhit Head 15m',
    'KL027':'Kalaloch 27m',
    'KL015':'Kalaloch 15m',
    'CE042':'Cape Elizabeth 42m',
    'CE015':'Cape Elizabeth 15m'  
}

mdates = pd.date_range('2011-01-01', periods=14, freq='YS')

# PLOTTING
plt.close('all')
fs=12
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,8))

# add mooring locations and labels to map 
axmap = plt.subplot2grid((5,3),(0,2),colspan=1,rowspan=5)
pfun.add_coast(axmap,color='grey')
pfun.dar(axmap)
axmap.axis([-125, -123, 47, 49])

axmap.text(-124.64,48.3,'Makah Bay',color='Black',weight='bold',alpha=0.8)
axmap.text(-124.67,48.13,'Cape Alava',color='Black',weight='bold',alpha=0.8)
axmap.text(-124.52,47.85,'Teahwhit Head',color='Black',weight='bold',alpha=0.8)
axmap.text(-124.35,47.57,'Kalaloch',color='Black',weight='bold',alpha=0.8)
axmap.text(-124.25,47.32,'Cape Elizabeth',color='Black',weight='bold',alpha=0.8)

axmap.set_xticks([-125,-124.5,-124,-123.5,-123])
axmap.set_xticklabels([-125,'',-124,'',-123])
axmap.set_yticks([47, 47.5, 48, 48.5, 49])
axmap.set_yticklabels([47, '', 48, '', 49])
axmap.set_xlabel('Longitude')
axmap.set_ylabel('Latitude')
axmap.set_title('OCNMS moorings')
                  
# plot data available for recent processed data 
sn_list = list(sn_name_dict.keys())
#sn_list = ['CE042']

ax1 = plt.subplot2grid((5,3),(0,0),colspan=1,rowspan=1)
ax2 = plt.subplot2grid((5,3),(0,1),colspan=1,rowspan=1)
ax3 = plt.subplot2grid((5,3),(1,0),colspan=1,rowspan=1)
ax4 = plt.subplot2grid((5,3),(1,1),colspan=1,rowspan=1)
ax5 = plt.subplot2grid((5,3),(2,0),colspan=1,rowspan=1)
ax6 = plt.subplot2grid((5,3),(2,1),colspan=1,rowspan=1)
ax7 = plt.subplot2grid((5,3),(3,0),colspan=1,rowspan=1)
ax8 = plt.subplot2grid((5,3),(3,1),colspan=1,rowspan=1)
ax9 = plt.subplot2grid((5,3),(4,0),colspan=1,rowspan=1)
ax10 = plt.subplot2grid((5,3),(4,1),colspan=1,rowspan=1)
 
# for Temp Sal marker 
marker_style = dict(color='navy', linestyle='none', marker='o',
                    markersize=10, markerfacecoloralt='tab:olive')
                  
ii = 1  
for sn in sn_list:
    print(sn)
    in_fn = data_dir / (sn + '_2011_2023_hourly.nc')
    
    ds = xr.open_dataset(in_fn, decode_times=True)
    IT = ds['IT'].values
    SAL = ds['SP'].values
    OXY = ds['DO (uM)'].values
    zii = ds.z.values
    
    # axes[0] is the sitemap 
    if ~np.isnan(IT).all():
        plt.gcf().axes[0].plot(ds.lon,ds.lat, fillstyle= 'full', **marker_style)
    if ~np.isnan(SAL).all():
        plt.gcf().axes[0].plot(ds.lon,ds.lat, fillstyle= 'top', **marker_style)
    if ~np.isnan(OXY).all():
        plt.gcf().axes[0].plot(ds.lon,ds.lat, marker = '+',color = 'tab:red',markersize=15)
                
    if ~np.isnan(OXY).all():
        o = OXY
        o[~np.isnan(OXY)]=1
        o=o*ds.z.values
        plt.gcf().axes[ii].plot(ds.time,o, marker = '.',color = 'tab:red',markersize=10)
    if ~np.isnan(IT).all():
        t = IT
        t[~np.isnan(IT)]=1
        t=t*ds.z.values
        plt.gcf().axes[ii].plot(ds.time,t, marker = '.',color = 'navy',markersize=3)    
    if ~np.isnan(SAL).all():
        s = SAL
        s[~np.isnan(SAL)]=1
        s=s*ds.z.values
        plt.gcf().axes[ii].plot(ds.time,s, marker = '.',color = 'tab:olive',markersize=3)
    
    plt.gcf().axes[ii].set_yticks(zii, minor=True)   
    if ii==7: 
        plt.gcf().axes[ii].set_ylim([-27,0])
        plt.gcf().axes[ii].set_yticks([0, -10, -20, -27], minor=False)
        plt.gcf().axes[ii].set_ylabel('Z est. [m]')
    elif ii%2: 
        plt.gcf().axes[ii].set_ylim([-42,0])
        plt.gcf().axes[ii].set_yticks([0, -10, -20, -30, -42], minor=False)
        plt.gcf().axes[ii].set_ylabel('Z est. [m]')
    else: 
        plt.gcf().axes[ii].set_ylim([-15,0])
        plt.gcf().axes[ii].set_yticks([0, -5, -10, -15], minor=False)
        
    plt.gcf().axes[ii].yaxis.grid(True, which='minor')
    plt.gcf().axes[ii].yaxis.grid(color = 'lightsteelblue', which='minor')
    plt.gcf().axes[ii].yaxis.grid(True, which='major')
    plt.gcf().axes[ii].yaxis.grid(color = 'Black', which='major')
          
    plt.gcf().axes[ii].set_xlim([mdates[0],mdates[-1]])   
    plt.gcf().axes[ii].set_xticks(mdates, minor=False)
    plt.gcf().axes[ii].xaxis.grid(True, which='major')   

    plt.gcf().axes[ii].tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False)
    
    #for label in plt.gcf().axes[ii].yaxis.get_ticklabels()[::2]:
    #    label.set_visible(False)
    
    plt.gcf().axes[ii].set_title(sn)
        
    #del IT, SAL, OXY 
    ii = ii + 1

for ii in [9, 10]:
    plt.gcf().axes[ii].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    plt.gcf().axes[ii].set_xticklabels(['J2011',' ','J2013',' ', 'J2015',' ',
                                        'J2017',' ','J2019',' ','J2021',' ','J2023',' '])
for ii in [1, 2]:  
    plt.gcf().axes[ii].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    plt.gcf().axes[ii].set_xticklabels([' ','J2012',' ','J2014', ' ','J2016',' ',
                                            'J2018',' ','J2020',' ','J2022',' ','J2024'])

plt.gcf().tight_layout()

axmap.plot(-124.7,48.81, fillstyle= 'full', **marker_style)
axmap.text(-124.7,48.81, "   temperature", horizontalalignment='left',verticalalignment='center', size=fs, color = 'navy',weight='bold')

axmap.plot(-124.64,48.755, fillstyle= 'top', **marker_style)
axmap.text(-124.64,48.755, "   salinity + temperature", horizontalalignment='left',verticalalignment='center', size=fs, color = 'tab:olive',weight='bold')

axmap.plot(-124.58,48.7, marker = '+',color = 'tab:red',markersize=10)
axmap.text(-124.58,48.7, "  dissolved oxygen", horizontalalignment='left',verticalalignment='center', size=fs, color = 'tab:red',weight='bold')

#fig_nm = data_dir / ('ocnms_2011_2023.png')
fig_nm = '/Users/katehewett/Documents/LKH_output/OCNMS_processed/ocnms_2011_2023.png'
plt.gcf().savefig(fig_nm)

#plt.show()

