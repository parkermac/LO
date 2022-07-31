"""
First code to compare observed to modeled casts.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

source = 'dfo'
year = 2019
year_str = str(year)
cbot = 'casts' # casts or bottles

# full DataFrame of observations
df = pd.read_pickle(Ldir['LOo'] / 'obs' / source / (cbot + '_' + str(year) + '.p'))

info_df = pd.read_pickle(Ldir['LOo'] / 'obs' / source / ('casts_' + str(year) + '.p'))

gtagex = 'cas6_v0_live'
mod_dir = Ldir['LOo'] / 'extract' / gtagex / 'cast' / (source + '_' + year_str)

plt.close('all')
pfun.start_plot(figsize=(20,8))

vn_dict = {'salt': 'salt (SA g kg-1)', 'temp': 'temp (CT degC)', 'oxygen': 'DO (uM)'}

ii = 0
for cid in info_df.index:
    mod_fn = mod_dir / (cbot + '_' + str(cid) + '.nc')
    if mod_fn.is_file() and (ii <= 20):
        ds = xr.open_dataset(mod_fn)
        cdf = df[df.cid==cid]
        fig = plt.figure()
        jj = 1
        for vn in vn_dict.keys():
            ax = fig.add_subplot(1,3,jj)
            ax.plot(ds[vn].values, ds.z_rho.values,'-r')
            ax.plot(cdf[vn_dict[vn]].to_numpy(), cdf['z'].to_numpy(),'-b')
            ax.set_title(vn)
            jj += 1
        ii += 1
        
plt.show()
pfun.end_plot()

