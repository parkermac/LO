"""
Make climatologies for point sources.
Discharge rate, temperature, and biogeochemisty variables.

Based on Ecology's timeseries, stored in LO_data/traps

To run, from ipython:
run make_climatology_pointsources.py
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates
import datetime


# helper function 
def monthly2daily(df):
    '''
    turn a monthly dataframe into daily data (for a year's worth of monthly data)
    '''
    # duplicate last row
    double_lr_df = pd.concat([df, df.iloc[-1:]], ignore_index=True)
    # picking arbitrary year to fill with daily data (but includes a leap year)
    start_date = datetime.date(2020, 1, 1)
    end_date = datetime.date(2021, 1, 1)
    dates = pd.date_range(start_date, end_date, freq='MS')
    # Replace month column with things that look like dates
    double_lr_df['Month'] = dates
    double_lr_df = double_lr_df.set_index('Month')
    # Change monthly to daily
    double_lr_daily_df = double_lr_df.resample('D').ffill()
    # delete last row (1/1 on the next year)
    daily_df = double_lr_daily_df[:-1]
    # make index start from 1 and go to 366
    daily_df.reset_index(inplace=True)
    return daily_df
    

# define year range to create climatologies
year0 = 1999
year1 = 2017

# location to save file
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'point_sources' /'Data_historical'

# file with all traps names and ID numbers
traps_info_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
# location of historical data to process
wwtp_dir = Ldir['data'] / 'traps' / 'point_sources'
all_his_fns = os.listdir(wwtp_dir)

# Get all wwtp names and wwtp IDs
traps_info_df = pd.read_excel(traps_info_fn,usecols='D,E,F')
# Only interested in wwtps
wwtp_all_df = traps_info_df.loc[traps_info_df['Inflow_Typ'] == 'Point Source']
#get river names
wwtp_names_df = wwtp_all_df['Name'].str.replace(' - 1', '')
# get river names and river ids
wwtpnames = wwtp_names_df.values
wwtpids = wwtp_all_df['ID'].values

# # just test Brightwater for now -------------------------------------------------
# wwtpnames = wwtpnames[46:47]
# wwtpids = wwtpids[46:47]

# # just test Birch Bay for now -------------------------------------------------
# wwtpnames = wwtpnames[26:27]
# wwtpids = wwtpids[26:27]

# # just test 5 WWTPs for now -------------------------------------------------
# wwtpnames = wwtpnames[12:17]
# wwtpids = wwtpids[12:17]

# initialize dataframes for all rivers
flow_clim_df = pd.DataFrame()
temp_clim_df = pd.DataFrame()
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
TIC_clim_df  = pd.DataFrame()
Talk_clim_df = pd.DataFrame()
DO_clim_df   = pd.DataFrame()

# variable names
vns = ['Flow(m3/s)','Temp(C)','NO3+NO2(mg/L)','NH4(mg/L)','DIC(mmol/m3)','Alk(mmol/m3)','DO(mg/L)']

# create one-year date range for plotting
yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')

# loop through all point sources
for i,wname in enumerate(wwtpnames):

    print('{}: {}'.format(i,wname))

    # get river index
    wID = wwtpids[i]
    
    wwtp_fn = ''

    # find Ecology's timeseries based on wwtp id
    for fn in all_his_fns:
        root, ext = os.path.splitext(fn)
        if root.startswith(str(wID)) and ext == '.xlsx':
            wwtp_fn = fn
    
    # Let user know if couldn't find timeseries for a given point source
    if wwtp_fn == '0':
        print('No history file found for {}'.format(wname))
    # Otherwise, read history file as df
    else:
        wwtp_fp = str(wwtp_dir)  + '/' + wwtp_fn
        wwtp_df = pd.read_excel(wwtp_fp, skiprows=[0])

    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    wwtp_df = wwtp_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)

    # replace all zeros with nans, so zeros don't bias data
    wwtp_df = wwtp_df.replace(0, np.nan)

    # calculate averages (compress 1999-2017 timeseries to single day, with an average for each day)
    wwtp_avgs_monthly_df = wwtp_df.groupby(['Month','Day']).mean().reset_index()
    # calculate standard deviation
    wwtp_sds_monthly_df = wwtp_df.groupby(['Month','Day']).std(ddof=0).reset_index()

    # convert monthly data to daily
    wwtp_avgs_df = monthly2daily(wwtp_avgs_monthly_df)
    wwtp_sds_df = monthly2daily(wwtp_sds_monthly_df)

    # replace all nans with zeros, so I'm no longer injecting nans
    wwtp_avgs_df = wwtp_avgs_df.replace(np.nan,0)
    wwtp_sds_df = wwtp_sds_df.replace(np.nan,0)

    # Set any negative TIC concentrations to zero
    wwtp_avgs_df['DIC(mmol/m3)'] = wwtp_avgs_df['DIC(mmol/m3)'].mask(wwtp_avgs_df['DIC(mmol/m3)'].lt(0),0)

    # Plot and save averages for each source
    fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
    ax = axes.ravel()
    for j,vn in enumerate(vns):
        i = j+1
        # label subplot
        ax[i].set_title(vn,fontsize=14)
        # Plot individual years
        for yr in range(1999,2017):
            wwtp_yr_monthly_df = wwtp_df.loc[wwtp_df['Year'] == yr]
            # convert monthly to daily
            wwtp_yr_df = monthly2daily(wwtp_yr_monthly_df)
            if yr == 2017:
                yrday_17 = pd.date_range(start ='1/1/2020', end ='8/02/2020', freq ='D') # don't have full 2017 dataset
                ax[i].plot(yrday_17,wwtp_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
            else:
                ax[i].plot(yrday,wwtp_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
        # Plot average
        ax[i].plot(yrday,wwtp_avgs_df[vn].values, label='average', color='black', linewidth=1.5)
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

    # Add data to climatology dataframes, and convert to units that LiveOcean expects
    flow_clim_df = pd.concat([flow_clim_df, pd.Series(wwtp_avgs_df['Flow(m3/s)'].values, name=wname)], axis = 1)     # [m3/s]
    temp_clim_df = pd.concat([temp_clim_df, pd.Series(wwtp_avgs_df['Temp(C)'].values, name=wname)], axis = 1)        # [C]
    NO3_clim_df  = pd.concat([NO3_clim_df, pd.Series(wwtp_avgs_df['NO3+NO2(mg/L)'] * 71.4, name=wname)], axis = 1)  # [mmol/m3]
    NH4_clim_df  = pd.concat([NH4_clim_df, pd.Series(wwtp_avgs_df['NH4(mg/L)'] * 71.4, name=wname)], axis = 1)      # [mmol/m3]
    TIC_clim_df  = pd.concat([TIC_clim_df, pd.Series(wwtp_avgs_df['DIC(mmol/m3)'], name=wname)], axis = 1)          # [mmol/m3]
    Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(wwtp_avgs_df['Alk(mmol/m3)'], name=wname)], axis = 1)          # [meq/m3]
    DO_clim_df   = pd.concat([DO_clim_df, pd.Series(wwtp_avgs_df['DO(mg/L)'] * 31.26, name=wname)], axis = 1)      # [mmol/m3]

# check for missing values:
if pd.isnull(flow_clim_df).sum().sum() != 0:
    print('Warning, there are missing flow values!')
if pd.isnull(temp_clim_df).sum().sum() != 0:
    print('Warning, there are missing temperature values!')
if pd.isnull(NO3_clim_df).sum().sum() != 0:
    print('Warning, there are missing nitrate values!')
if pd.isnull(NH4_clim_df).sum().sum() != 0:
    print('Warning, there are missing ammonium values!')
if pd.isnull(TIC_clim_df).sum().sum() != 0:
    print('Warning, there are missing TIC values!')
if pd.isnull(Talk_clim_df).sum().sum() != 0:
    print('Warning, there are missing alkalinity values!')
if pd.isnull(DO_clim_df).sum().sum() != 0:
    print('Warning, there are missing oxygen values!')

# save results
flow_clim_df.to_pickle(clim_dir / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p'))
temp_clim_df.to_pickle(clim_dir / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p'))
NO3_clim_df.to_pickle(clim_dir / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p'))
NH4_clim_df.to_pickle(clim_dir / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p'))
TIC_clim_df.to_pickle(clim_dir / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p'))
Talk_clim_df.to_pickle(clim_dir / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p'))
DO_clim_df.to_pickle(clim_dir / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p'))

# Calculate summary statistics for all wwtps
clim_avgs = pd.DataFrame()
clim_max = pd.DataFrame()
clim_min = pd.DataFrame()
clim_sds = pd.DataFrame()
# list climatology dfs
clim_df_list = [flow_clim_df,temp_clim_df,NO3_clim_df,NH4_clim_df,TIC_clim_df,Talk_clim_df,DO_clim_df]
for i,vn in enumerate(vns):
    scale = 1
    if vn == 'NO3+NO2(mg/L)' or vn == 'NH4(mg/L)':
        scale = 71.4
    elif vn == 'DO(mg/L)':
        scale = 31.26
    # average values of all point sources
    clim_avgs[vn] = clim_df_list[i].mean(axis=1)/scale
    # max climatology values
    clim_max[vn] = clim_df_list[i].max(axis=1)/scale
    # min climatology values
    clim_min[vn] = clim_df_list[i].min(axis=1)/scale
    # standard deviation of all point sources
    clim_sds[vn] = clim_df_list[i].std(axis=1)/scale

# Plot Summary Statistics
fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
ax = axes.ravel()
for j,vn in enumerate(vns):
    i = j+1
    # label subplot
    ax[i].set_title(vn,fontsize=14)
    # Plot average
    ax[i].plot(yrday,clim_avgs[vn].values, label='Average of all Sources', color='mediumpurple', linewidth=1.5)
    # Plot error shading
    upper_bound = [min(clim_avgs[vn].values[ii]+clim_sds[vn].values[ii],clim_max[vn].values[ii]) for ii in range(366)] # don't go higher than max value
    lower_bound = [max(clim_avgs[vn].values[ii]-clim_sds[vn].values[ii],clim_min[vn].values[ii]) for ii in range(366)] # don't go lower than min value
    ax[i].fill_between(yrday,upper_bound,lower_bound,label='One SD',color='mediumpurple',alpha=0.2,edgecolor='none')
    # Plot max
    ax[i].plot(yrday,clim_max[vn].values, label='Max Value', color='firebrick', linestyle='--', linewidth=1)
    # Plot min
    ax[i].plot(yrday,clim_min[vn].values, label='Min Value', color='cornflowerblue', linestyle='--', linewidth=1)
    # fontsize of tick labels
    ax[i].tick_params(axis='both', which='major', labelsize=12)
    ax[i].tick_params(axis='x', which='major', rotation=30)
    ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
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
figname = 'point_source_summary.png'
save_path = clim_dir / figname
fig.savefig(save_path)
# plt.close('all')
plt.show()
