"""
Plot as-run river time series.

"""
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart(gridname='cas6', tag='v3')

# load extraction (an xarray Dataset)
fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms' / 'extraction_2017.01.01_2021.12.31.nc'
x = xr.load_dataset(fn)

# get climatology
clm_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_historical' / 'CLIM_flow_1980_2020.p'
dfc = pd.read_pickle(clm_fn)

# add the climatology, for practice
x['transport_clim'] = 0*x.transport
x['yearday'] = (('time'), x.time.to_index().dayofyear.to_numpy())
ydvec = x.yearday.values

# add the climatology to the xarray dataset (maybe use groupby instead?)
for rn in list(x.riv.values):
    if rn in dfc.columns:
        this_riv = dfc[rn] # a Series
        this_riv_clim = 0 * ydvec
        for ii in range(1,367):
            this_riv_clim[ydvec==ii] = this_riv[ii]
        x.transport_clim.loc[:,rn] = this_riv_clim
    else:
        print('Missing ' + rn)
        
# plotting
plt.close('all')
pfun.start_plot()
fig = plt.figure()

# for plotting time series we are better off using pandas
df = pd.DataFrame(index=x.time.values)
ii = 1
for rn in ['fraser', 'columbia', 'skagit', 'deschutes']:
    ax = fig.add_subplot(2,2,ii)
    df.loc[:,'Q'] = x.transport.sel(riv=rn).values
    df.loc[:,'Qclim'] = x.transport_clim.sel(riv=rn).values
    if ii == 1:
        leg = True
    else:
        leg = False
    df.plot(ax=ax, grid=True, legend=leg)
    ax.set_ylim(bottom=0)
    ax.set_xlim(x.time.values[0], x.time.values[-1])
    ax.text(.05,.9,rn.title(), transform=ax.transAxes)
    if ii in [1,2]:
        ax.set_xticklabels([])
    ii += 1

plt.show()
pfun.end_plot()
