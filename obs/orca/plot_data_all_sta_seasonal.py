"""
Code to plot the processed ORCA mooring data, combining all stations
on a single plot, and showing seasonal profiles of one property per plot.
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

sn_list = ['CI','PW','NB','DB','HP','TW']
axnum_dict = dict(zip(sn_list,[1,2,4,5,7,8]))

axlim_dict = {'SA':(24,31),'CT':(8,18),'DO (uM)':(0,400),'SIG0':(18,24)}

if testing:
    vn_list = ['SA']
else:
    vn_list = ['SA','CT','DO (uM)','SIG0']
    
plt.close('all')
pfun.start_plot(figsize=(16,8))

for vn in vn_list:
    print( vn)
    fig = plt.figure()

    jj = 1
    for sn in sn_list:
        print(' - '+sn)
        out_fn = out_dir / (sn + '_daily.nc')
        ds = xr.open_dataset(out_fn)

        z = ds.z.values
        t = ds.time.to_index()
    
        fld = np.transpose(ds[vn].values) # pack as (z,time) for plotting
        
        # map
        if jj == 1:
            axm = fig.add_subplot(133)
            pfun.add_coast(axm)
            axm.axis([-123.5,-122,47,48.3])
            pfun.dar(axm)
            axm.set_title('ORCA Station Location')
            axm.set_xlabel('Longitude')
            axm.set_ylabel('Latitude')
        axm.plot(ds.attrs['lon'],ds.attrs['lat'],'or',ms=10,alpha=.3)
        axm.text(ds.attrs['lon'],ds.attrs['lat'],sn+': '+ds.attrs['Station Name'],
            bbox=pfun.bbox,rotation=20)
    
        # time-mean profiles, by season
        ax = fig.add_subplot(3,3,axnum_dict[sn])
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
            if jj == 1:
                ax.text(.95,.35 - ii*.1,obfun.season_name_dict[ii],color=c_dict[ii],
                    fontweight='bold',ha='right',transform=ax.transAxes)
        ax.text(.05,.15,sn,color='k',
            fontweight='bold',transform=ax.transAxes)
        ax.set_ylim(-120,0)
        ax.set_xlim(axlim_dict[vn][0], axlim_dict[vn][1])
        ax.grid(True)
        # ax.set_xlim(xmin=0)
        if jj in [5,6]:
            ax.set_xlabel('Seasonal Mean ' + ds[vn].attrs['long_name']+' ['+ds[vn].attrs['units']+']')
        if jj in [1,3,5]:
            ax.set_ylabel('Z [m]')

        fig.tight_layout()
        
        jj += 1
    
    if testing:
        plt.show()
    else:
        out_name = 'All_sta_seasonal_' + vn.replace(' ','_') + '.png'
        plt.savefig(plot_out_dir / out_name)
        plt.close()
