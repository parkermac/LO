"""
This is code for making the time series expected by the "superplot" series
in plotting/roms_plots.py.

Run on mac in ipython:
run make_superplot_forcing -gtx cas6_v3_lo8b -y 2019

"""
import sys
import xarray as xr
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt

from lo_tools import Lfun, zrfun, zfun

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
# select time period and frequency
parser.add_argument('-y', '--year', type=int) # e.g. 2019

# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided (other)
argsd = args.__dict__
for a in ['gtagex', 'year']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

gtagex = Ldir['gtagex']
year = str(Ldir['year'])
d_str = year + '.01.01_' + year + '.12.31'
        
# define where to put the output
outdir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'superplot'
Lfun.make_dir(outdir)

# Rivers
riv_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms' / ('extraction_' + d_str + '.nc')
riv_ds = xr.open_dataset(riv_fn)
# get DataArray of transport for only selected rivers
flow_da = riv_ds.transport.sel(riv=['columbia', 'fraser', 'skagit'])

# Mooring
moor_fn = Ldir['LOo'] / 'extract' / gtagex / 'moor' / ('superplot_' + d_str + '.nc')
moor_ds = xr.open_dataset(moor_fn)
eta_rms = np.sqrt(zfun.lowpass(moor_ds.zeta.values**2, f='godin'))[1::24]
svstr_lp = zfun.filt_AB8d(moor_ds['svstr'].values)[1::24]

# Combined
comb_df = pd.DataFrame(index=pd.date_range(start='1/1/'+year, end='12/31/'+year))
comb_df.loc[:,'RMS Tide Height (m)'] = eta_rms
comb_df.loc[:,'8-day NS Wind Stress (Pa)'] = svstr_lp
comb_df.loc[:,'Columbia R. Flow (1000 m3/s)'] = flow_da.sel(riv='columbia')/1000
comb_df.loc[:,'Fraser R. Flow (1000 m3/s)'] = flow_da.sel(riv='fraser')/1000
comb_df.loc[:,'Skagit R. Flow (1000 m3/s)'] = flow_da.sel(riv='skagit')/1000

# save for later use
out_fn = outdir / ('forcing_' + gtagex + '_' + year + '.p')
print('Saving ' + str(out_fn))
comb_df.to_pickle(out_fn)

# plotting
plt.close('all')
comb_df.plot(subplots=True, xlim=(comb_df.index[0], comb_df.index[-1]), grid=True)
plt.show()

