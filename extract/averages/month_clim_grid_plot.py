"""
Code to explore plotting of monthly means and their anomaly from climatology.

This version is intended to put all the fields or anomalies on a single grid.
The goal is to make it easier to see patterns.

To get the full suite of plots the main things to change are -pt (mean or anomaly)
and -vn (temp or oxygen).

It assumes that you want to plot all years and all months that are available
in the files in the LO_roms/[gtagex]/averages and climatologies folders that
you have. It does this by processing the directory contents.

"""

import xarray as xr
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from datetime import datetime, timedelta

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

parser = argparse.ArgumentParser()

# command line arguments
parser.add_argument('-gtx', '--gtagex', default='cas7_t0_x4b', type=str)
parser.add_argument('-ro', '--roms_out_num', default=0, type=int)
# parser.add_argument('-0', '--ds0', default='2014.01.01', type=str)
# parser.add_argument('-1', '--ds1', default = '2023.12.31', type=str)
parser.add_argument('-pt', '--plot_type', default='mean', type=str) # mean or anomaly
parser.add_argument('-vn', default='temp', type=str) # temp or oxygen
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

# do things with the arguments
args = parser.parse_args()
argsd = args.__dict__
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# dt0 = datetime.strptime(Ldir['ds0'], Lfun.ds_fmt)
# dt1 = datetime.strptime(Ldir['ds1'], Lfun.ds_fmt)
# dti = pd.date_range(dt0, dt1, freq='ME', inclusive='both')

dir1 = Ldir['roms_out'] / Ldir['gtagex'] / 'averages'
dir2 = Ldir['roms_out'] / Ldir['gtagex'] / 'climatologies'
fn1_list = list(dir1.glob('*.nc'))
fn1_list.sort()

# year list repeats each year 12 times
year_list = [item.name.split('.')[0].split('_')[-2] for item in fn1_list]

# month list repeats 01-12 as many times as there are years
month_list = [item.name.split('.')[0].split('_')[-1] for item in fn1_list]

fn2_list = list(dir2.glob('*.nc'))
fn2_list.sort()

# fn1_list = []
# fn2_list = []
# year_list = []
# month_list = []
# for dt in dti:
#     ym_str = dt.strftime('%Y_%m')
#     mo_str = ('000' + str(dt.month))[-2:]
#     fn1_list.append(dir1 / ('monthly_mean_' + ym_str + '.nc'))
#     fn2_list.append(dir2 / ('monthly_clim_' + mo_str + '.nc'))
#     year_list.append(dt.year)
#     month_list.append(dt.month)

# find out how many years we have, to help format the plots
year_list_unique = list(set(year_list))
year_list_unique.sort() # in place
nyears = len(year_list_unique)

# plotting choices

# color lims for field and difference from climatology
c_dict = {'temp': (4,20), 'oxygen':(0,400)}
d_dict = {'temp': (-3,3), 'oxygen':(-40,40)}

m_dict = dict(zip(list(range(1,13)),list('JFMAMJJASOND')))

# PLOT CODE
vn = Ldir['vn']
if vn == 'temp':
    tstr = 'Surface Temp. (deg C)'
    slev = -1 # surface plot
elif vn == 'oxygen':
    tstr = 'Bottom DO (uM)'
    slev = 0 # bottom plot

pt = Ldir['plot_type']
if pt == 'mean':
    pstr = 'Monthly Mean'
    vmin, vmax = c_dict[vn]
elif pt == 'anomaly':
    pstr = 'Monthly Anomaly'
    vmin, vmax = d_dict[vn]

# plotting

plt.close('all')
if Ldir['testing']:
    # small plot for laptop
    pfun.start_plot(figsize=(6,8))
else:
    pfun.start_plot(figsize=(8,12))
fig = plt.figure()
gs1 = gridspec.GridSpec(nyears,12, figure=fig)
gs1.update(wspace=0.025, hspace=0.015) # set the spacing between axes. 

ii = 0
for fn1 in fn1_list:

    if ii == 0:
        ds1 = xr.open_dataset(fn1)
        plon, plat = pfun.get_plon_plat(ds1['lon_rho'].values, ds1['lat_rho'].values)
        f1 = ds1[vn][0,slev,:,:].values
        ds1.close()
        if pt == 'anomaly':
            fn2 = fn2_list[ii]
            ds2 = xr.open_dataset(fn2)
            f2 = ds2[vn][0,slev,:,:].values
            ds2.close()

    # if Ldir['testing'] and (ii > 0):
    #     pass
    # elif (not Ldir['testing']) and (ii > 0):
    else:
        ds1 = xr.open_dataset(fn1)
        f1 = ds1[vn][0,slev,:,:].values
        ds1.close()
        if pt == 'anomaly':
            fn2 = fn2_list[ii]
            ds2 = xr.open_dataset(fn2)
            f2 = ds2[vn][0,slev,:,:].values
            ds2.close()

    ax1 = plt.subplot(gs1[ii])
    if pt == 'mean':
        cmap='rainbow'
        cs1 = ax1.pcolormesh(plon,plat,f1,cmap=cmap, vmin=vmin, vmax=vmax)
    elif pt == 'anomaly':
        cmap='RdYlBu_r'
        cs1 = ax1.pcolormesh(plon,plat,f1-f2,cmap=cmap, vmin=vmin, vmax=vmax)
    pfun.dar(ax1)
    # ax1.axis('off')
    if ii in range(12):
        ax1.set_title(m_dict[month_list[ii]])
    if ii in range(0,nyears*12,12):
        ax1.set_ylabel(year_list[ii])
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_xticks([])
    ax1.set_yticks([])

    ii += 1

fig.suptitle('%s: %s: (Range = %d, %d)' % (tstr, pstr, vmin, vmax))

out_dir = Ldir['LOo'] / 'plots'
out_fn = out_dir / (vn + '_' + pt + '.png')
print(out_fn)

if Ldir['testing']:
    plt.show()
else:
    Lfun.make_dir(out_dir)
    fig.savefig(out_fn)
