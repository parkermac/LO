"""
Code to test the module loadDFO.py
"""
from datetime import datetime, timedelta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from lo_tools import plotting_functions as pfun
import sys

import loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# this is where the data files are stored
basedir=str(Ldir['data'] / 'obs' / 'dfo')

# make some choices about time and space range
year = 2017
datelims=[datetime(year,1,1),datetime(year,12,31)]
latlims=[48.9,49.5]
lonlims=[-124.2,-123.2]

mo_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
    'August', 'September', 'October', 'November', 'December']
mo_dict = dict(zip(range(1,13), mo_list))

plt.close('all')

# make some colors
colors = plt.cm.hsv(np.linspace(0,1,12))
colors = np.roll(colors, 6, axis=0)

if False:
    # CTD Data
    c_df = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO_CTD.sqlite',
            datelims=datelims,latlims=latlims,lonlims=lonlims)
    # plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    for icast in c_df['Station'].unique():
        ax.plot(c_df.loc[c_df.Station==icast,['SA']],
                   c_df.loc[c_df.Station==icast,['Z']],'-')
    ax.set_ylabel('Z (m)')

if False:
    # CTD Data, locations only
    c_df = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO_CTD.sqlite',
        datelims=(),latlims=(),lonlims=(), xyt_only=True)
    # plot
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111)
    c_df.plot(x='Lon',y='Lat',style='.', ax=ax, legend=False)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis([-132,-122,46,53])
    
    # Print number of casts in each year.
    # RESULTS: There are only casts in these year:
    # 2013: 1077
    # 2014: 995
    # 2015: 1216
    # 2016: 1307
    # 2017: 1128
    # 2018: 593
    # 2019: 1188
    # 2020: 55
    for year in range(1900,datetime.now().year+1):
        df = c_df[c_df.dtUTC.dt.year == year]
        nsta = len(df.Station.unique())
        if nsta>0:
            print('%d: %d' % (year, nsta))

if True:
    # Bottle Data
    # b_df = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
    #     datelims=datelims,latlims=latlims,lonlims=lonlims)
    
    year = 2017
    datelims=[datetime(year,1,1),datetime(year,12,31)]
    latlims=[48.8,50.0]
    lonlims=[-125,-122.5]
    
    b_df = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
        datelims=datelims, lonlims=lonlims, latlims=latlims)
    if b_df is None:
        print('no data!')
        sys.exit()
    # plot
    fig = plt.figure(figsize=(18,10))
    vn_list = ['SA','CT','N','DO', 'Chl', 'Si']
    ii = 1
    for vn in vn_list:
        if vn == 'SA':
            units = 'g/kg'
        elif vn == 'CT':
            units = 'degC'
        else:
            a = b_df[vn+'_units'].unique()
            b = [item for item in a if item != None]
            if len(b) > 1:
                print('Warning: extra units for ' + vn)
                print(b)
            else:
                units = b[0]
        ax = fig.add_subplot(2,3,ii)
        for mo in range(1,13):
            df = b_df[b_df['dtUTC'].dt.month==mo]
            if len(df) > 0:
                for icast in df['Station'].unique():
                    df[df.Station==icast].plot(x=vn,y='Z', ax=ax, style='*', color=colors[mo-1], legend=False)
            else:
                print('no data for month = ' + str(mo))
        # add month lables
        if ii == 1:
            for mo in range(1,13):
                ax.text(.05, .95 - .06*mo, mo_dict[mo], color=colors[mo-1], fontweight='bold',
                    transform=ax.transAxes)
        ii += 1
        ax.set_xlabel(vn + ' [' + units + ']')
        ax.set_ylabel('Z [m]')
        
        
if False:
    # Bottle Data, locations only
    b_df = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
        datelims=(),latlims=(),lonlims=(), xyt_only=True)
    # plot
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111)
    b_df.plot(x='Lon',y='Lat',style='.', ax=ax, legend=False)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis([-132,-122,46,53])
    
    # Print number of casts in each year.
    # RESULTS: the bottle data goes 1930-2019).  Pretty solid after 2000, but
    # 1932, -53, and -54 were pretty impressive as well. There are stations
    # in Puget Sound but they look like they have a longitude offset, often
    # sitting on land.
    for year in range(1900,datetime.now().year+1):
        df = b_df[b_df.dtUTC.dt.year == year]
        nsta = len(df.Station.unique())
        if nsta>0:
            print('%d: %d' % (year, nsta))
    


plt.show()


