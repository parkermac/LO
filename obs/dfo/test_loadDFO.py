"""
Code to test the module loadDFO.py
"""
from datetime import datetime, timedelta
from pathlib import Path
import matplotlib.pyplot as plt

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


c_df = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO_CTD.sqlite',
        datelims=datelims)
        
b_df = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
    datelims=datelims,latlims=latlims,lonlims=lonlims)

# plotting
plt.close('all')

fig,ax=plt.subplots(1,2,figsize=(10,10))
fig.subplots_adjust(wspace=.5)

for icast in c_df['CTDStationTBLID'].unique():
    ax[0].plot(c_df.loc[c_df.CTDStationTBLID==icast,['SA']],
               c_df.loc[c_df.CTDStationTBLID==icast,['Z']],'k-')
ax[0].set_xlabel('Reference Salinity (g kg$^{-1}$)')
ax[0].set_ylabel('Depth (m)')
ax[0].set_title('March 2016, \n Central Strait Of Georgia')

for icast in b_df['BOTStationTBLID'].unique():
    ax[1].plot(b_df.loc[b_df.BOTStationTBLID==icast,['N']],
               b_df.loc[b_df.BOTStationTBLID==icast,['Z']],'k-*')
#ax[1].set_ylim((350,0))
ax[1].set_xlabel('Nitrate+Nitrite (mmol N m$^{-3}$)')
ax[1].set_ylabel('Depth (m)')
ax[1].set_title('2016, \n Central Strait Of Georgia')

plt.show()


