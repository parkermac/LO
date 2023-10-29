"""
Code to plot the processed ORCA mooring data, combining all stations
on a single plot, and showing time series of three properties.

For the LO/obs README.
"""

import pandas as pd
import xarray as xr
import numpy as np
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

plot_out_dir = Ldir['LO'] / 'obs' / 'readme_plots' # place for README figures
Lfun.make_dir(plot_out_dir)
plot_out_fn = plot_out_dir / 'orca.png'

sn_list = ['TW','HP','CI','PW','NB','DB'] # move TW first because it is noisy
sn_c_dict = {'CI':'g','PW':'b','NB':'r','DB':'salmon','HP':'orange','TW':'gold'}

if testing:
    sn_list = ['CI']

axlim_dict = {'SA':(24,31),'CT':(8,18),'DO (uM)':(0,400)}

vn_list = ['SA','CT','DO (uM)']
axnum_dict = {'SA':1,'CT':3,'DO (uM)':5}

plt.close('all')
pfun.start_plot(figsize=(18,8))
fig = plt.figure()

ss = 1
ax_dict = dict()
for sn in sn_list:
    print(sn)
    out_fn = out_dir / (sn + '_daily.nc')
    
    ds = xr.open_dataset(out_fn)

    z = ds.z.values
    t = ds.time.to_index()

    vv = 1
    for vn in vn_list:
        print(' - '+vn)
        
        fld = np.transpose(ds[vn].values) # pack as (z,time) for plotting
        
        # map
        if (vv == 1) and (ss == 1):
            axm = fig.add_subplot(133)
            pfun.add_coast(axm)
            axm.axis([-123.5,-122,47,48.3])
            pfun.dar(axm)
            axm.set_title('ORCA Station Location')
            axm.set_xlabel('Longitude')
            axm.set_ylabel('Latitude')
        if vv == 1:
            axm.plot(ds.attrs['lon'],ds.attrs['lat'],'o',color=sn_c_dict[sn],ms=10)
            axm.text(.05+ds.attrs['lon'],ds.attrs['lat'],sn+': '+ds.attrs['Station Name'],
                bbox=pfun.bbox,rotation=20)
    
        # time series
        if ss == 1:
            ax_dict[vv] = plt.subplot2grid((3,3), (vv-1,0), colspan=2)
        ax = ax_dict[vv]
            
        # mask out depths with poor time coverage makes jumpy averages
        nans = np.isnan(fld).sum(axis=1)
        nans_norm = nans/np.min(nans)
        zmask = nans_norm > 1.25
        fld[zmask,:] = np.nan
        
        # average over a depth range
        zdiv = -10
        mask = z >= zdiv
        this_fld = fld[mask,:]
        fld_mean = np.nanmean(this_fld, axis=0)
        # monthly means to make plots clearer
        ser = pd.Series(index=t, data=fld_mean)
        serm = ser.resample('M').mean()
        st = serm.index
        sv = serm.to_numpy()
        # plot line for this time series
        ax.plot(st,sv,'-',color=sn_c_dict[sn],lw=3)
        
        if ss == 1:
            ax.text(.05,.15,vn,color='k',
                fontweight='bold',transform=ax.transAxes,bbox=pfun.bbox)
            ax.set_xlim(pd.Timestamp(2005,1,1), pd.Timestamp(2021,12,31))
            # ax.set_ylim(axlim_dict[vn][0], axlim_dict[vn][1])
            ax.grid(True)
            if vv == 3:
                ax.set_xlabel('Time')
            if vv in [1,2]:
                ax.set_xticklabels([])
            if vv == 1:
                ax.set_title('Monthly means above Z = %d [m]' % (zdiv))

        
        vv += 1 # vn
        
    ss += 1 # station
    
fig.tight_layout()
plt.show()

if not testing:
    plt.savefig(plot_out_fn)
