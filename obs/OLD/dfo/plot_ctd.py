"""
Code to look at the output of ctd processing.

Note, it appears that individual casts are packed bottom top to bottom.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

year = 2017
aa = [-125.2,-122.5, 48.5,50.0]
x0, x1, y0, y1 = aa

df = pd.read_pickle(Ldir['LOo'] / 'obs' / 'dfo' / 'ctd' / (str(year) + '.p'))
# ['sta', 'lon', 'lat', 'time', 'z', 'SA', 'CT', 'DO (uM)', 'Fluor']

# keep only data in a box
df = df[(df['lon']>x0) & (df['lon']<x1) & (df['lat']>y0) & (df['lat']<y1)]

# trim ocean sections as well with a diagonal line
df = df[df['lat'] > (-.8/1.2)*(df['lon']+125.2) + 49.4]

plt.close('all')
pfun.start_plot()

fig = plt.figure(figsize=(22,12))

vn_list = ['SA', 'CT', 'DO (uM)', 'Fluor']
ax_list = [1,2,5,6]
ax_dict = dict(zip(vn_list,ax_list))

# map axis
axm = fig.add_subplot(122)

for vn in vn_list:
    ii = ax_dict[vn]
    ax = fig.add_subplot(2,4,ii)
    for mo in range(1,13):
        dfm = df[df['time'].dt.month==mo]
        if len(dfm) > 0:
            dfm.plot(x=vn,y='z', ax=ax, style='.', color=pfun.month_color_dict[mo], legend=False)
            if vn == 'SA':
                dfm.plot(x='lon',y='lat', ax=axm, style='o',
                    color=pfun.month_color_dict[mo], legend=False,
                    markersize=26-2*mo)
    # add month labels
    if ii == 1:
        for mo in range(1,13):
            ax.text(.05, .95 - .065*mo, pfun.month_name_dict[mo],
                color=pfun.month_color_dict[mo], fontweight='bold',
                transform=ax.transAxes)
    ax.set_xlabel('')
    ax.text(.95,.05,vn,transform=ax.transAxes,ha='right',fontweight='bold')
    if ii in [1,5]:
        ax.set_ylabel('Z [m]')
    
pfun.add_coast(axm)
pfun.dar(axm)
axm.axis(aa)
axm.set_title('DFO SoG CTD Data ' + str(year))
axm.set_xlabel('Longitude')
axm.set_ylabel('Latitude')

plt.show()
pfun.end_plot()
