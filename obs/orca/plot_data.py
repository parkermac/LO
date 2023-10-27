"""
Code to plot the processed ORCA mooring data, mainly as a
quick check to see if the fields look reasonable.
"""

import pandas as pd
import xarray as xr
import numpy as np
import gsw
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings('ignore') # skip warning about all-nan layers

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

source = 'orca'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle
out_dir = Ldir['LOo'] / 'obs' / source / otype

if testing:
    sn_list = ['CI']
else:
    sn_list = ['CI','PW','NB','DB','HP','TW']
    
plt.close('all')
pfun.start_plot(figsize=(20,10))

for sn in sn_list:
    out_fn = out_dir / (sn + '_daily.nc')
    ds = xr.open_dataset(out_fn)
    
    fig = plt.figure()
    
    z = ds.z.values
    t = ds.time.to_index()
    
    # make extended axes for plotting
    Zp = z[1:] + np.diff(z)/2
    Zm = z[:-1] - np.diff(z)/2
    Z = np.concatenate((Zm[:2],Zp))
    Tp = t + timedelta(days=0.5)
    Tm = t - timedelta(days=0.5)
    T = Tm.union(Tp)
    
    vn = 'DO (uM)'
    fld = np.transpose(ds[vn].values)
    
    # time series
    ax = fig.add_subplot(221)
    cs = ax.pcolormesh(T,Z,fld)
    fig.colorbar(cs, ax=ax)
    ax.set_title('ORCA Mooring %s: %s' % (ds.attrs['Station Name'],vn))
    ax.set_xlabel('Time')
    ax.set_ylabel('Z [m]')
    ax.set_ylim(Z[0],0)
    
    # map
    ax = fig.add_subplot(122)
    ax.plot(ds.attrs['lon'],ds.attrs['lat'],'or',ms=10)
    pfun.add_coast(ax)
    ax.axis([-123.5,-122,47,48.3])
    pfun.dar(ax)
    ax.set_title('Station Location')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    # time-mean profiles, by season
    ax = fig.add_subplot(223)
    c_dict = {0:'b',1:'g',2:'r',3:'orange'}
    daylist = [0, 92,183,275,367]
    for ii in range(len(daylist)-1):
        tmask = (t.dayofyear>=daylist[ii]) & (t.dayofyear<daylist[ii+1])
        this_fld = fld[:,tmask]
        # mask out depths where poor time coverage makes jumpy averages
        nans = np.isnan(this_fld).sum(axis=1)
        nans_norm = nans/np.min(nans)
        fld_mean = np.nanmean(this_fld,axis=1)
        fld_mean[nans_norm > 1.25] = np.nan
        
        ax.plot(fld_mean,z,'-',color=c_dict[ii],lw=3)
    ax.set_ylim(Z[0],0)
    ax.grid(True)
    # ax.set_xlim(xmin=0)
    ax.set_xlabel('Time Mean ' + vn)
    ax.set_ylabel('Z [m]')
    
    
    
plt.show()
