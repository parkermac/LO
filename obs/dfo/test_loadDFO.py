"""
Code to test the module loadDFO.py
"""
from datetime import datetime, timedelta
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

import loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# this is where teh data files are stored
basedir=str(Ldir['data'] / 'obs' / 'dfo')

# make some choices about time and space range
datelims=[datetime(2016,1,1),datetime(2016,12,31)]
latlims=[48.9,49.5]
lonlims=[-124.2,-123.2]

plt.close('all')

# make some colors
colors = plt.cm.rainbow(np.linspace(0,1,12))

if False:
    # CTD Data
    c_df = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO_CTD.sqlite',
            datelims=datelims)
    # plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    for icast in c_df['CTDStationTBLID'].unique():
        ax.plot(c_df.loc[c_df.CTDStationTBLID==icast,['SA']],
                   c_df.loc[c_df.CTDStationTBLID==icast,['Z']],'-')
    ax.set_xlabel('Reference Salinity (g kg$^{-1}$)')
    ax.set_ylabel('Z (m)')
    ax.set_title('CTD 2016, \n Strait Of Georgia')

if True:
    # Bottle Data
    b_df = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
        datelims=datelims,latlims=latlims,lonlims=lonlims)
    # plot
    fig = plt.figure(figsize=(14,10))
    vn_list = ['SA','CT','N','DO']
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
        ax = fig.add_subplot(2,2,ii)
        for mo in range(1,13):
            df = b_df[b_df['dtUTC'].dt.month==mo]
            for icast in df['Station'].unique():
                df[df.Station==icast].plot(x=vn,y='Z', ax=ax, style='-*', color=colors[mo-1], legend=False)
        ii += 1
        ax.set_xlabel(vn + ' [' + units + ']')
        ax.set_ylabel('Z [m]')

plt.show()


