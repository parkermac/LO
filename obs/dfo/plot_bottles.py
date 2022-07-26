"""
Code to look at the output of process_bottles.py.
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

df = pd.read_pickle(Ldir['LOo'] / 'obs' / 'dfo' / ('bottles_' + str(year) + '.p'))
# ['sta', 'lon', 'lat', 'time', 'z', 'salt (SA g kg-1)', 'temp (CT degC)', 'DO (uM)', 'NO3 (uM)']

# keep only data in a box
df = df[(df['lon']>x0) & (df['lon']<x1) & (df['lat']>y0) & (df['lat']<y1)]

# trim ocean sections as well with a diagonal line
df = df[df['lat'] > (-.8/1.2)*(df['lon']+125.2) + 49.4]

months = range(1,13) # a list of 1 to 12
month_name_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
    'August', 'September', 'October', 'November', 'December']
month_name_dict = dict(zip(range(1,13), month_name_list))
month_color_dict = dict(zip(months,
    ['mediumblue', 'royalblue', 'cadetblue', 'aquamarine',
    'lightgreen', 'greenyellow', 'gold', 'orange',
    'lightsalmon', 'mediumorchid', 'slateblue', 'purple']))


plt.close('all')
pfun.start_plot()

fig = plt.figure(figsize=(22,10))

vn_list = ['salt (SA g kg-1)', 'temp (CT degC)', 'DO (uM)', 'NO3 (uM)']
ax_list = [1,2,5,6]
ax_dict = dict(zip(vn_list,ax_list))

# map axis
axm = fig.add_subplot(122)

for vn in vn_list:
    ii = ax_dict[vn]
    ax = fig.add_subplot(2,4,ii)
    for mo in months:
        dfm = df[df['time'].dt.month==mo]
        if len(dfm) > 0:
                dfm.plot(x=vn,y='z', ax=ax, style='.', color=month_color_dict[mo], legend=False)
                dfm.plot(x='lon',y='lat', ax=axm, style='o', color=month_color_dict[mo], legend=False,
                markersize=26-2*mo)#, markeredgecolor='k')
    # add month lables
    if ii == 1:
        for mo in range(1,13):
            ax.text(.05, .95 - .06*mo, month_name_dict[mo],
                color=month_color_dict[mo], fontweight='bold',
                transform=ax.transAxes)
    ax.set_xlabel(vn)
    if ii in [1,5]:
        ax.set_ylabel('Z [m]')
    
pfun.add_coast(axm)
pfun.dar(axm)
axm.axis(aa)
axm.set_title('DFO SoG Bottle Data ' + str(year))
axm.set_xlabel('Longitude')
axm.set_ylabel('Latitude')


plt.show()
pfun.end_plot()
