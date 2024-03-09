"""
Make climatologies for point sources.
Discharge rate, temperature, and biogeochemisty variables.

Based on Ecology's timeseries, using data stored in 
LO_data/[trapsD##]/all_point_source_data.nc
(To change the Ecology data version, modify traps_data_ver.csv)

        To run, from ipython:
        run make_climatology_pointsources.py

To create individual climatology figures, run from ipython with:
run make_climatology_pointsources.py -test True

Figures saved in:
LO_output/pre/trapsP##/point_sources/[ctag]/Data_historical/climatology_plots

Note that running with -test True adds
several minutes to run time. (~ 3 min)
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import datetime
import matplotlib.dates as mdates
import datetime
import traps_helper
import os
from pathlib import Path
    

#################################################################################
#                     Get data and set up dataframes                            #
#################################################################################

# read arguments
parser = argparse.ArgumentParser()
# -test True will output plots
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# add ctag
parser.add_argument('-ctag', type=str, default='lo_base')
args = parser.parse_args()
ctag = args.ctag

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent
with open(this_dir / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# get traps pre-processing version (folder name)
path = os.getcwd()
trapsP_vers = os.path.basename(path)

# location to save file
clim_dir = Ldir['LOo'] / 'pre' / trapsP_vers / 'point_sources' / ctag / 'Data_historical'
Lfun.make_dir(clim_dir)

# get flow and loading data
wwtp_fn = Ldir['data'] / trapsD / 'all_point_source_data.nc'
ecology_data_ds = xr.open_dataset(wwtp_fn)

# get wwtp names and wwtp ids
wwtpnames = ecology_data_ds['name'].values

# # just test a few WWTPs for now 
# wwtpnames = wwtpnames[28:30]

# initialize dataframes for all wwtps
DO_clim_df   = pd.DataFrame()
flow_clim_df = pd.DataFrame()
temp_clim_df = pd.DataFrame()
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
TIC_clim_df  = pd.DataFrame()
Talk_clim_df = pd.DataFrame()

# variable names
vns = ['DO(mmol/m3)','Flow(m3/s)','Temp(C)','NO3(mmol/m3)',
       'NH4(mmol/m3)','TIC(mmol/m3)','Talk(meq/m3)']
clims = [DO_clim_df, flow_clim_df, temp_clim_df, NO3_clim_df,
         NH4_clim_df, TIC_clim_df, Talk_clim_df]
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']

# create one-year date range for plotting
yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')

#################################################################################
#                          Calculate climatologies                              #
#################################################################################

print('Calculating point source climatologies...')

# loop through all point sources
for i,wname in enumerate(wwtpnames):

    # turn dataset information for this wwtp into a dataframe
    wwtp_df, wwtp_avgs_df, wwtp_sds_df = traps_helper.ds_to_avgdf(wname,ecology_data_ds)

#################################################################################
#                            Create climatologies                               #
#################################################################################

    # Add data to climatology dataframes
    for j,vn in enumerate(vns):
        clims[j] = pd.concat([clims[j],pd.Series(wwtp_avgs_df[vn].values, name=wname)],
                         axis = 1)

# Calculate summary statistics for all wwtps
clim_avgs = pd.DataFrame()
clim_max = pd.DataFrame()
clim_min = pd.DataFrame()
clim_sds = pd.DataFrame()

for i,vn in enumerate(vns):
    # average values of all point sources
    clim_avgs[vn] = clims[i].mean(axis=1)
    # max climatology values
    clim_max[vn] = clims[i].max(axis=1)
    # min climatology values
    clim_min[vn] = clims[i].min(axis=1)
    # standard deviation of all point sources
    clim_sds[vn] = clims[i].std(axis=1)

#################################################################################
#                           Plot summary statistics                             #
#################################################################################

# Plot Summary Statistics
fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
ax = axes.ravel()
for j,vn in enumerate(vns):
    i = j+1

    # convert DO from mmol/m3 to mg/L for plotting
    if vn == 'DO(mmol/m3)':
        scale = 1/31.26 
        var = 'DO(mg/L)'
    else:
        scale = 1
        var = vn

    # get average values, standard deviation, min and max
    avgs =     clim_avgs[vn].values * scale
    sds =      clim_sds[vn].values  * scale
    maxvals =  clim_max[vn].values  * scale
    minvals =  clim_min[vn].values  * scale

    # label subplot
    ax[i].text(0.05,0.85,letters[j]+' '+var,transform=ax[i].transAxes,fontsize=14)
    # Plot average
    ax[i].plot(yrday,avgs, label='Average of all Sources',
               color='mediumpurple', linewidth=1.5)
    # Plot error shading
    upper_bound = [min(avgs[ii]+sds[ii],maxvals[ii]) for ii in range(366)] # don't go higher than max value
    lower_bound = [max(avgs[ii]-sds[ii],minvals[ii]) for ii in range(366)] # don't go lower than min value
    ax[i].fill_between(yrday,upper_bound,lower_bound,label='One SD',
                       color='mediumpurple',alpha=0.2,edgecolor='none')
    # Plot max
    ax[i].plot(yrday,maxvals, label='Max Value',
               color='firebrick', linestyle='--', linewidth=1)
    # Plot min
    ax[i].plot(yrday,minvals, label='Min Value',
               color='cornflowerblue', linestyle='--', linewidth=1)
    # fontsize of tick labels
    ax[i].tick_params(axis='both', which='major', labelsize=12)
    ax[i].tick_params(axis='x', which='major', rotation=30)
    ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
    if i < 7:
        ax[i].set_ylim([0,1.3*max(maxvals)])
    # create legend
    if i ==7:
        ax[i].set_ylim([0,1.3*max(clim_max['TIC(mmol/m3)'].values)])
        handles, labels = ax[7].get_legend_handles_labels()
        ax[0].legend(handles, labels, loc='center', ncol = 2,fontsize=14)
        ax[0].axis('off')
    # Define the date format
    if i >= 6:
        date_form = mdates.DateFormatter("%b")
        ax[i].xaxis.set_major_formatter(date_form)
# plot title is name of source
plt.suptitle('Point Source Climatology Summary (n={})'.format(len(wwtpnames)),fontsize=18)
# Save figure
figname = 'point_source_summary.png'
save_path = clim_dir / figname
fig.savefig(save_path)

print('Point source climatologies done\n')

#################################################################################
#              Plot all individual point source climatologies                   #
#################################################################################

if args.testing == True:
    print('Plotting...')

    # generate directory to save files
    fig_dir = clim_dir / 'climatology_plots'
    Lfun.make_dir(fig_dir)

    for i,wname in enumerate(wwtpnames):

        print('{}/{}: {}'.format(i+1,len(wwtpnames),wname))

        # turn dataset information for this wwtp into a dataframe
        wwtp_df, wwtp_avgs_df, wwtp_sds_df = traps_helper.ds_to_avgdf(wname,ecology_data_ds)

        # Plot climatologies for each source
        fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
        ax = axes.ravel()
        for j,vn in enumerate(vns):

            # convert DO from mmol/m3 to mg/L for plotting
            if vn == 'DO(mmol/m3)':
                scale = 1/31.26 
                var = 'DO(mg/L)'
            else:
                scale = 1
                var = vn

            i = j+1
            # label subplot
            ax[i].set_title(var,fontsize=14)
            # Plot individual years
            for yr in range(1999,2017):
                wwtp_yr_df = wwtp_df.loc[wwtp_df['year'] == yr]
                values_to_plot = wwtp_yr_df[vn].values*scale
                values_to_plot = values_to_plot.tolist()
                # skip leap years
                if np.mod(yr,4) != 0:
                    # pad Feb 29th with nan
                    values_to_plot = values_to_plot[0:60] + [np.nan] + values_to_plot[60::]
                if yr == 2017:
                    yrday_17 = pd.date_range(start ='1/1/2020',
                                            end ='8/02/2020', freq ='D') # don't have full 2017 dataset
                    ax[i].plot(yrday_17,values_to_plot,alpha=0.5, label=yr, linewidth=1)
                else:
                    ax[i].plot(yrday,values_to_plot,alpha=0.5, label=yr, linewidth=1)
            # Plot average
            ax[i].plot(yrday,wwtp_avgs_df[vn].values*scale, label='average', color='black', linewidth=1.5)
            # fontsize of tick labels
            ax[i].tick_params(axis='both', which='major', labelsize=12)
            ax[i].tick_params(axis='x', which='major', rotation=30)
            ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
            # create legend
            if i ==7:
                handles, labels = ax[7].get_legend_handles_labels()
                ax[0].legend(handles, labels, loc='center', ncol = 4,fontsize=14)
                ax[0].axis('off')
            # Define the date format
            if i >= 6:
                date_form = mdates.DateFormatter("%b")
                ax[i].xaxis.set_major_formatter(date_form)
        # plot title is name of source
        plt.suptitle(wname,fontsize=18)
        # Save figure
        figname = wname + '.png'
        save_path = clim_dir / 'climatology_plots' / figname
        fig.savefig(save_path)
        plt.close('all')

    print('Done')
      
#################################################################################
#                             Save climatologies                                #
#################################################################################

pickle_names = ['DO', 'flow', 'temp', 'NO3', 'NH4', 'TIC', 'Talk']

# check for missing values:
for i,clim in enumerate(clims):
    if pd.isnull(clim).sum().sum() != 0:
        print('Warning, there are missing '+pickle_names[i]+' values!')

# save results
for i,clim in enumerate(clims):
    clim.to_pickle(clim_dir / ('CLIM_' + pickle_names[i] + '.p'))
    # # test printing
    # print(vns[i]+'=================================================')
    # print(clim)