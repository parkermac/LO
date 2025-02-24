"""
General-purpose code to look at the output of the ctd processing for a given
source and year.

Example command to run in ipython:
run plot_ctd_bottle -source nceiSalish -otype bottle -year 2017

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-source', type=str, default='ecology_nc') # e.g. dfo1
parser.add_argument('-otype', type=str, default='bottle') # observation type, e.g. ctd, bottle
parser.add_argument('-year', type=int, default = 2017) # e.g. 2019
parser.add_argument('-small', type=Lfun.boolean_string, default = False) # True for laptop size

args = parser.parse_args()

if args.source=='nceiSalish':
    aa = [-125.5, -122, 47, 50]
elif args.source in ['ecology','ecology_nc']:
    aa = [-125, -122, 46, 50]
else:
    aa = [-130, -122, 42, 52]
x0, x1, y0, y1 = aa

df = pd.read_pickle(Ldir['LOo'] / 'obs' / args.source / args.otype / (str(args.year) + '.p'))

# keep only data in a box
df = df[(df['lon']>x0) & (df['lon']<x1) & (df['lat']>y0) & (df['lat']<y1)]
zbot = df.z.min()

plt.close('all')

if args.small:
    pfun.start_plot(fs=10)
    fig = plt.figure(figsize=(12,8)) # laptop size
else:
    pfun.start_plot()
    fig = plt.figure(figsize=(22,13))

if args.otype == 'bottle':
    vn_list = ['SA', 'CT', 'DO (uM)', 'Chl (mg m-3)', 'NO3 (uM)', 'NH4 (uM)', 'DIC (uM)', 'TA (uM)']
    ax_list = [1,2,4,5,7,8,10,11]
    nrows = 4
    ncols = 3
    markersize = 5
elif args.otype == 'ctd':
    vn_list = ['SA', 'CT', 'DO (uM)', 'Chl (mg m-3)']
    ax_list = [1,2,4,5]
    nrows = 2
    ncols = 3
    markersize = 1
else:
    print('Unexpected otype: ' + args.otype)
    sys.exit()
ax_dict = dict(zip(vn_list,ax_list))

# map axis
axm = fig.add_subplot(1,ncols,ncols)

n_dict = {}

for vn in vn_list:
    if vn in df.columns:
        ii = ax_dict[vn]
        ax = fig.add_subplot(nrows,ncols,ii)
        for mo in range(1,13):
            dfm = df[df['time'].dt.month==mo]
            if vn == vn_list[0]:
                n_dict[mo] = len(dfm) # keep track of number of samples in each month
            if len(dfm) > 0:
                dfm.plot(x=vn,y='z', ax=ax, style='.',
                    color=pfun.month_color_dict[mo], legend=False,
                    markersize=markersize)
                if vn == vn_list[0]:
                    dfm.plot(x='lon',y='lat', ax=axm, style='o',
                        color=pfun.month_color_dict[mo], legend=False,
                        markersize=26-2*mo)
        # add month labels and number of samples
        if vn == vn_list[0]:
            for mo in range(1,13):
                ax.text(.05, .97 - .077*mo, '%d %s' % (n_dict[mo], pfun.month_name_dict[mo]),
                    color=pfun.month_color_dict[mo],# fontweight='bold',
                    transform=ax.transAxes, ha = 'left')
        ax.set_xlabel('')
        ax.set_ylim(bottom=zbot-20)
        # ax.set_xlim(left=0)
        ax.text(.95,.05,vn,transform=ax.transAxes,ha='right',fontweight='bold')
        if ii in [1,5,9]:
            ax.set_ylabel('Z [m]')
    
pfun.add_coast(axm)
pfun.dar(axm)
axm.axis(aa)
axm.set_title('%s %s Data for %d' % (args.source, args.otype.upper(), args.year))
axm.set_xlabel('Longitude')
axm.set_ylabel('Latitude')

fig.tight_layout()
plt.show()
pfun.end_plot()
