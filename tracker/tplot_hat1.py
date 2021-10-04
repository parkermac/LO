"""
Plot results of a particle tracking experiment for the man overboard
near Hat Island from November 2020.

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
sys.path.append(os.path.abspath('../plotting'))
import pfun
import seawater as sw

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
import pytz

Ldir = Lfun.Lstart()

# Set an experiment to plot from.
indir0 = Ldir['LOo'] + 'tracks/'
outdir = indir0 + 'hat_plots_2021_05/'
Lfun.make_dir(outdir, clean=True)

exp_list = ['hat_surf_sh4','hat_surf_sh5','hat_surf_sh6',
            'hat_surf_wind1_sh4','hat_surf_wind1_sh5','hat_surf_wind1_sh6',
            'hat_surf_wind2_sh4','hat_surf_wind2_sh5','hat_surf_wind2_sh6']

for exp in exp_list:
    indir = exp + '/'
    rel = 'release_2020.11.18.nc'

    # get Datasets
    fn = indir0 + indir + rel
    dsr = nc.Dataset(fn)

    NT, NP = dsr['lon'].shape

    # get a list of datetimes
    ot_vec = dsr['ot'][:]
    dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

    def get_dt_local(dt, tzl='US/Pacific'):
        tz_utc = pytz.timezone('UTC')
        tz_local = pytz.timezone(tzl)
        dt_utc = dt.replace(tzinfo=tz_utc)
        dt_local = dt_utc.astimezone(tz_local)
        return dt_local
    
    # tracks (packed [time, particle])
    lon = dsr['lon'][:]
    lat = dsr['lat'][:]
    h = dsr['h'][:]
    zeta = dsr['zeta'][:,0]
    x0 = lon[0,:]
    y0 = lat[0,:]

    dtl_list = [get_dt_local(dt) for dt in dt_list]
    h_ser = pd.Series(index=dtl_list, data=zeta)

    # PLOTTING
    plt.close('all')
    fs = 16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(16,10))

    # MAP
    ax = fig.add_subplot(121)
    aa = [-122.7, -122.1, 47.8, 48.4]
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-122.6, -122.4, -122.2])
    ax.set_yticks([47.8, 48, 48.2, 48.4])
    ax.text(.05, .05, exp, transform=ax.transAxes, weight='bold', c='gray')

    c_list = ['r','orange','gold','y','limegreen','steelblue','b']

    # TIME SERIES
    ax2 = fig.add_subplot(122)
    ax2.set_xlabel('Date and Time [PST]')
    ax2.set_ylabel('Sea Surface Height [m]')
    
    # add the tracks (packed [time, particle]) colored by time segments
    i0 = 0
    step=6
    i1 = i0 + step
    cc = 0
    while i1 <= NT:
        ax.plot(lon[i0:i1+1,:], lat[i0:i1+1,:], '.-', c=c_list[cc], lw=1, alpha=.3)
        ax.plot(lon[0,:].mean(), lat[0,:].mean(), '*k', markersize=20)
        h_ser.iloc[i0:i1+1].plot(ax=ax2, linestyle='-', marker='o', c=c_list[cc], lw=3, grid=True)
        i0 += step
        i1 = i0 + step
        cc += 1
    ax.plot(lon[0,:].mean(), lat[0,:].mean(), '*k', markersize=20)

    fig.tight_layout()
    fig.savefig(outdir + exp + '.png')
    #plt.show()
    plt.rcdefaults()

    dsr.close()

