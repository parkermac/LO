"""
Make climatologies for river flow and temperature. Also saves plots.

This code shows how powerful pandas is for this kind of task.
Really it is just one line to make a climatology (the groupby call)
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-ctag', type=str, default='lo_base')
args = parser.parse_args()
ctag = args.ctag

# location of historical data to plot
riv_dir0 = Ldir['LOo'] / 'pre' / 'river1' / ctag
riv_dir = riv_dir0 / 'Data_historical'
plot_dir = riv_dir0 / 'Data_historical_plots'
Lfun.make_dir(plot_dir)

# flow
flow_df = pd.read_pickle(riv_dir / ('ALL_flow.p'))
# make the climatology (one line!)
flow_clim_df = flow_df.groupby(flow_df.index.dayofyear).mean()
# check for missing values:
if pd.isnull(flow_clim_df).sum().sum() != 0:
    print('Warning, there are missing flow values!')
# save results
flow_clim_df.to_pickle(riv_dir / ('CLIM_flow.p'))

# temperature
temp_df = pd.read_pickle(riv_dir / ('ALL_temperature.p'))
temp_clim_df = temp_df.groupby(temp_df.index.dayofyear).mean()
# drop temperature rivers with more than 50 missing yeardays
temp_clim_df = temp_clim_df.loc[:, pd.isnull(temp_clim_df).sum() < 50]
# fill missing temperature values
temp_clim_df = temp_clim_df.interpolate()
# make yearday 366 equal to yearday 365 (leap year is poorly sampled)
temp_clim_df.loc[366,:] = temp_clim_df.loc[365,:]
if pd.isnull(temp_clim_df).sum().sum() != 0:
    print('Warning, there are missing temperature values!')
temp_clim_df.to_pickle(riv_dir / ('CLIM_temp.p'))

# Plotting
plt.close('all')

fig = plt.figure(figsize=(16,10))
rn_split = np.array_split(flow_clim_df.columns, 9)
for ii in range(1,10):
    ax = fig.add_subplot(3,3,ii)
    flow_clim_df[rn_split[ii-1]].plot(ax=ax)
    ax.set_xlim(0,366)
    ax.set_ylim(bottom=0)
    if ii >= 7:
        ax.set_xlabel('Yearday')
    if ii in [1, 4, 7]:
        ax.set_ylabel(r'Flow [$m^{3}s^{-1}$]')
fig.savefig(plot_dir / ('CLIM_flow_plot.png'))

fig = plt.figure(figsize=(16,10))
rn_split = np.array_split(temp_clim_df.columns, 9)
for ii in range(1,10):
    ax = fig.add_subplot(3,3,ii)
    temp_clim_df[rn_split[ii-1]].plot(ax=ax)
    ax.set_xlim(0,366)
    ax.set_ylim(0, 25)
    if ii >= 7:
        ax.set_xlabel('Yearday')
    if ii in [1, 4, 7]:
        ax.set_ylabel(r'Temperature [$^{\circ}C$]')
fig.savefig(plot_dir / ('CLIM_temp_plot.png'))
    
plt.show()