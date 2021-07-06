"""
Plot as-run river time series.

"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta

import Lfun
import zrfun
import zfun

Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms' / 'extraction_2018.01.01_2018.12.31.p'
df = pd.read_pickle(fn)

# get climatology
dfi = df.index

cols = df.columns

clm_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_historical' / 'CLIM_flow_1980_2020.p'
df_clim = pd.read_pickle(clm_fn)

rind = df.index
dt0 = rind[0]
dt1 = rind[-1]

# plotting
plt.close('all')

for ff in range(5): # range(5)
    
    fig = plt.figure(figsize=(12,10))
    
    for ii in range(9):
        rr = ii + 9*ff
        riv_name = df.iloc[:,rr].name
        ax = fig.add_subplot(3,3,ii+1)
        df.iloc[:,rr].plot(ax=ax, color='purple')
        df_clim.iloc[:,rr].plot(ax=ax, color='orange')
        ax.set_xlim(dt0, dt1)
        ax.text(.05, .9, riv_name, transform=ax.transAxes, fontweight='bold')
        for year in range(2018,2021):
            ax.axvline(x=datetime(year,1,1), c='k')
        ax.set_xticks([], minor=True)
        if ii+1 <=6:
            ax.set_xticklabels([])
    fig.suptitle(fnr)
        
plt.show()
