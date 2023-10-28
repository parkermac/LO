"""
Code to plot the processed ORCA mooring data, mainly as a check to see
if the fields look reasonable.
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
from lo_tools import obs_functions as obfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

source = 'orca'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle
out_dir = Ldir['LOo'] / 'obs' / source / otype

plot_out_dir = out_dir / 'figs' # general archive of figures
Lfun.make_dir(plot_out_dir)

if testing:
    sn_list = ['CI']
    vn_list = ['DO (uM)']
else:
    sn_list = ['CI','PW','NB','DB','HP','TW']
    vn_list = ['SA','CT','DO (uM)','SIG0']
    
plt.close('all')
pfun.start_plot(figsize=(16,8))

for sn in sn_list:
    print(sn)
    out_fn = out_dir / (sn + '_daily.nc')
    ds = xr.open_dataset(out_fn)
    
    for vn in vn_list:
        print(' - ' + vn)
        fig = plt.figure()
    
        z = ds.z.values
        t = ds.time.to_index()
    
        # make extended axes for pcolormesh
        Zp = z[1:] + np.diff(z)/2
        Zm = z[:-1] - np.diff(z)/2
        Z = np.concatenate((Zm[:2],Zp))
        Tp = t + timedelta(days=0.5)
        Tm = t - timedelta(days=0.5)
        T = Tm.union(Tp) # cute
    
        fld = np.transpose(ds[vn].values) # pack as (z,time) for plotting
    
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
        # get season info dicts
        c_dict = obfun.season_c_dict
        daylist = obfun.season_daylist
        for ii in range(len(daylist)-1): # loop over seasons 0-3
            tmask = (t.dayofyear>=daylist[ii]) & (t.dayofyear<daylist[ii+1])
            this_fld = fld[:,tmask]
            # mask out depths where poor time coverage makes jumpy averages
            nans = np.isnan(this_fld).sum(axis=1)
            nans_norm = nans/np.min(nans)
            fld_mean = np.nanmean(this_fld,axis=1)
            fld_mean[nans_norm > 1.25] = np.nan
            # plot line for this season
            ax.plot(fld_mean,z,'-',color=c_dict[ii],lw=3)
            # name of season
            ax.text(.95,.35 - ii*.1,obfun.season_name_dict[ii],color=c_dict[ii],
                fontweight='bold',ha='right',transform=ax.transAxes)
        ax.text(.05,.15,'Seasonal Mean Profiles',color='k',
            fontweight='bold',transform=ax.transAxes)
        ax.text(.05,.05,ds[vn].attrs['long_name']+' ['+ds[vn].attrs['units']+']',
            color='k',fontweight='bold',transform=ax.transAxes)
        ax.set_ylim(Z[0],0)
        ax.grid(True)
        # ax.set_xlim(xmin=0)
        ax.set_xlabel('Time Mean ' + vn)
        ax.set_ylabel('Z [m]')

        fig.tight_layout()
    
        if testing:
            plt.show()
        else:
            out_name = 'Single_property_' + sn + '_' + vn.replace(' ','_') + '.png'
            plt.savefig(plot_out_dir / out_name)
            plt.close()
