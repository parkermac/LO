"""
Plot as-run river time series.

"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import plotting_functions as pfun

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart(gridname='cas6', tag='v3')

# load extraction (an xarray Dataset)
fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms' / 'extraction_2018.01.01_2018.12.31.nc'
x = xr.load_dataset(fn)

# get climatology
clm_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_historical' / 'CLIM_flow_1980_2020.p'
dfc = pd.read_pickle(clm_fn)

# add the climatology, for practice
x = x.assign(transport_clim = 0*x.transport)


# add the climatology to the xarray dataset
NT = len(x.time) # needed because the climatology has 366 days
for rn in list(x.riv.values):
    if rn in dfc.columns:
        x.transport_clim.loc[dict(riv=rn)] = dfc[rn][:NT].to_numpy()
    else:
        print('Missing ' + rn)
        
# plotting
plt.close('all')
pfun.start_plot()
fig = plt.figure()

# for time series we are better off using pandas
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
