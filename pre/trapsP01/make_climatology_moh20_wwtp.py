"""
Make climatologies for Mohamedali et al. (2020) WWTPs.
Discharge rate, temperature, and biogeochemisty variables.

Based on Mohamedali et al. (2020) timeseries, using data stored in 
LO_data/trapsD01/processed_data/wwtp_data_mohamedali_etal_2020.nc

        To run, from ipython:
        run make_climatology_moh20_wwtp.py

To create individual climatology figures, run from ipython with:
run make_climatology_moh20_wwtp.py -test True

Figures saved in:
LO_output/pre/trapsP01/wwtps/[ctag]/Data_historical/climatology_plots

Note that running with -test True adds
several minutes to run time. (~ 1 min)
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

plt.close('all')

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
clim_dir = Ldir['LOo'] / 'pre' / trapsP_vers / 'moh20_wwtps' / ctag / 'Data_historical'
Lfun.make_dir(clim_dir)

# get flow and loading data
moh20_wwtp_fn = Ldir['data'] / trapsD / 'processed_data'/ 'wwtp_data_mohamedali_etal_2020.nc'
moh20_wwtp_data_ds = xr.open_dataset(moh20_wwtp_fn)

# get wwtp names and wwtp ids
wwtpnames = moh20_wwtp_data_ds['name'].values

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

# initialize variable name and climatology df dictionary
vn_clim_dict = {'Flow(m3/s)': flow_clim_df,
                'NO3(mmol/m3)': NO3_clim_df,
                'NH4(mmol/m3)': NH4_clim_df,
                'TIC(mmol/m3)': TIC_clim_df,
                'Talk(meq/m3)': Talk_clim_df,
                'DO(mmol/m3)': DO_clim_df,
                'Temp(C)': temp_clim_df}

# initialize pickle name dictionary
vn_pickle_dict = {'Flow(m3/s)': 'flow',
                'NO3(mmol/m3)': 'NO3',
                'NH4(mmol/m3)': 'NH4',
                'TIC(mmol/m3)': 'TIC',
                'Talk(meq/m3)': 'Talk',
                'DO(mmol/m3)': 'DO',
                'Temp(C)': 'temp'}

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
    wwtp_df, wwtp_avgs_df, wwtp_sds_df = traps_helper.ds_to_avgdf(wname,moh20_wwtp_data_ds)

#################################################################################
#                            Create climatologies                               #
#################################################################################

    # Add WWTP's data as new column to climatology dataframe
    for vn, clim_df in vn_clim_dict.items():
        new_col = pd.Series(wwtp_avgs_df[vn].values, name=wname)
        vn_clim_dict[vn] = pd.concat([clim_df, new_col], axis=1)


# Calculate summary statistics for all wwtps
clim_avgs = pd.DataFrame()
clim_max = pd.DataFrame()
clim_min = pd.DataFrame()
clim_sds = pd.DataFrame()

for vn in vn_clim_dict.keys():
    # average values of all point sources
    clim_avgs[vn] = vn_clim_dict[vn].mean(axis=1)
    # max climatology values
    clim_max[vn] = vn_clim_dict[vn].max(axis=1)
    # min climatology values
    clim_min[vn] = vn_clim_dict[vn].min(axis=1)
    # standard deviation of all point sources
    clim_sds[vn] = vn_clim_dict[vn].std(axis=1)
    

#################################################################################
#                           Plot summary statistics                             #
#################################################################################

# Plot Summary Statistics
fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
ax = axes.ravel()
for j,vn in enumerate(vn_clim_dict.keys()):
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
    ax[i].set_ylim([0,1.3*max(maxvals)])
    # create legend
    if i ==7:
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
figname = 'moh20_wwtp_summary.png'
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

    for j,wname in enumerate(wwtpnames):

        print('{}/{}: {}'.format(j+1,len(wwtpnames),wname))

        # turn dataset information for this wwtp into a dataframe
        wwtp_df, wwtp_avgs_df, wwtp_sds_df = traps_helper.ds_to_avgdf(wname,moh20_wwtp_data_ds)

        # Plot climatologies for each source
        fig, axes = plt.subplots(7,1, figsize=(10, 9), sharex=True)
        ax = axes.ravel()

        # loop through variables
        for i,vn in enumerate(vn_clim_dict.keys()):

            # format axis
            ax[i].set_facecolor('#EEEEEE')
            ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
            for border in ['top','right','bottom','left']:
                ax[i].spines[border].set_visible(False)

            # convert DO from mmol/m3 to mg/L for plotting
            if vn == 'DO(mmol/m3)':
                scale = 1/31.26 
                var = 'DO(mg/L)'
            else:
                scale = 1
                var = vn

            # label subplot
            ax[i].text(0.02,0.8,var,fontsize=12,fontweight='bold',
                       transform=ax[i].transAxes)

            # Plot whole time series
            values_to_plot = wwtp_df[vn].values*scale # scale variable properly
            time = wwtp_df['Date'].values
            ax[i].plot(time,values_to_plot, color='deeppink', label='Raw data',
                       linewidth=2, alpha=0.5)

            # Plot average
            leapyear_clim = wwtp_avgs_df[vn].values*scale
            # remove leap day to get non-leap year climatology (365 days)
            feb29_index = 59
            nonleapyear_clim = np.delete(leapyear_clim, feb29_index)
            # generate climatology for the entire period
            # correctly alternate leap and non-leap years to create artificial 2005-2020 time series
            leap = leapyear_clim   # renamed for convenience
            non = nonleapyear_clim # renamed for convenience
            clim_1999to2017 = np.concatenate((
            #   1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017
                non ,leap,non ,non ,non ,leap,non ,non ,non ,leap,non ,non ,non ,leap,non ,non ,non ,leap,non[0:212]))
            # plot
            ax[i].plot(time,clim_1999to2017, label='climatology',
                        color='black', linewidth=1)
            
            # set y-axis
            ax[i].set_ylim([0,1.3*np.nanmax(values_to_plot)])
            
            # add legend
            if i == 0:
                ax[i].legend(loc='upper right',ncol=2,fontsize=12)

            # fontsize of tick labels
            ax[i].tick_params(axis='both', which='major', labelsize=12)
            ax[i].tick_params(axis='x', which='major', rotation=30)
            ax[i].set_xlim([datetime.date(1999, 1, 1), datetime.date(2017, 12, 31)])

            # Define the date format
            ax[1].xaxis.set_major_locator(mdates.YearLocator())
            ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

        # plot title is name of source
        plt.suptitle(wname,fontsize=14)
        plt.subplots_adjust(top=0.95)

        # Save figure
        figname = wname + '.png'
        save_path = clim_dir / 'climatology_plots' / figname
        fig.savefig(save_path)
        plt.close('all')

    print('Done')
      
#################################################################################
#                             Save climatologies                                #
#################################################################################

# check for missing values:
for vn,clim in vn_clim_dict.items():
    if pd.isnull(clim).sum().sum() != 0:
        print('Warning, there are missing '+vn+' values!')

# save results
for vn,pickle_name in vn_pickle_dict.items():
    clim = vn_clim_dict[vn]
    clim.to_pickle(clim_dir / ('CLIM_' + pickle_name + '.p'))