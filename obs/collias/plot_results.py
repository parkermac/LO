"""
Code to plot the Collias data, after running process_data.py.

The idea is to check visually on the data.

"""

import pandas as pd
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing
small = args.small

# BOTTLE
source = 'collias'
otype = 'bottle'
in_dir = Ldir['data'] / 'obs' / source

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype

# plot output location
if not testing:
    plot_out_dir = Ldir['LOo'] / 'obs' / source / (otype+'_check_plots')
    Lfun.make_dir(plot_out_dir)

# default lists
year_list = range(1932,1976)

if testing:
    print_info = True
    year_list = [1932]
    # sta_list = ['LCH551']
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.options.display.width = 0 # auto-detect full display width

# Plotting
plt.close('all')
if small:
    pfun.start_plot(figsize=(12,8))
else:
    pfun.start_plot(figsize=(20,13))

for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')

    if out_fn.is_file():
        df = pd.read_pickle(out_fn)
    else:
        print('- no data')
        continue
        
    sta_list = list(df.name.unique())
    
    # checking for bad values
    if testing:
        df1 = df.loc[df['NO2 (uM)']>30,:]
        if len(df1)>0:
            print('** Bad data found in year '+ys)

    fig = plt.figure()

    # map
    ax = fig.add_subplot(331)
    df.plot(x='lon', y='lat',style='.b', legend=False, ax=ax)
    if testing:
        df1.plot(x='lon', y='lat',style='.r', legend=False, ax=ax)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-125, -122, 47, 49])
    ax.set_title('Collias Dataset ' + ys)
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    # time series
    ax = plt.subplot2grid((3,3), (0,1), colspan=2)
    df.plot(x='time',y='z',style='.b', legend=False, ax=ax)
    if testing:
        df1.plot(x='time',y='z',style='.r', legend=False, ax=ax)
    ax.set_xlabel('Time')
    ax.set_ylabel('Z [m]')
    ax.grid(True)
    # tracers
    ii = 4
    vn_list = ['CT', 'SA', 'DO (uM)','NO3 (uM)', 'NO2 (uM)', 'SiO4 (uM)']
    for vn in vn_list:
        if vn in df.columns:
            ax = fig.add_subplot(3,3,ii)
            df.plot(x=vn, y='z',style='.b', legend=False, ax=ax, grid=True)
            if testing:
                df1.plot(x=vn, y='z',style='.r', legend=False, ax=ax, grid=True)
            if ii in [4,7]:
                ax.set_ylabel('Z [m]')
            ax.text(.95,.1,vn,fontweight='bold',transform=ax.transAxes,ha='right')
        else:
            pass
        ii += 1
        
    fig.tight_layout()
    if testing:
        plt.show()
    else:
        plt.savefig(plot_out_dir / ('check_plot_'+ys+'.png'))
        plt.close()
        
pfun.end_plot()
