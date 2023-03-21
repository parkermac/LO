"""
Make climatologies for tiny rivers.
Discharge rate, temperature, and biogeochemisty variables.

Based on Ecology's timeseries, stored in LO_data/traps

To run, from ipython:
run make_climatology_tinyrivs.py
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates

# define year range to create climatologies
year0 = 1999
year1 = 2017

# location to save files
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'tiny_rivers' /'Data_historical'

# file with all traps names and ID numbers
traps_info_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
# location of historical data to process
riv_dir = Ldir['data'] / 'traps' / 'nonpoint_sources'
all_his_fns = os.listdir(riv_dir)

# Get all river names and river IDs
traps_info_df = pd.read_excel(traps_info_fn,usecols='D,E,F')
# Only interested in rivers
riv_all_df = traps_info_df.loc[traps_info_df['Inflow_Typ'] == 'River']
# remove double counted river names (some are listed in two rows)
riv_singles_df = riv_all_df.loc[traps_info_df['Name'].str.contains('- 2') == False]
# rename rivers that have a ' - 1' at the end
riv_withLOrepeats_df = riv_singles_df.replace(' - 1', '', regex=True)
# read overlapping rivers
repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
repeatrivs_df = pd.read_excel(repeatrivs_fn)
SSM_repeats = repeatrivs_df['SSM_rname'].values
# remove repeat river names from list of river names
riv_noLOrepeats_df = riv_withLOrepeats_df[~riv_withLOrepeats_df['Name'].isin(SSM_repeats)]
# get river names and river ids
rivnames = riv_noLOrepeats_df['Name'].values
rivids = riv_noLOrepeats_df['ID'].values

# # just Tsitsika River for now -------------------------------------------------
# rivnames = ['Tsitika']
# rivids = [448]

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

# list of rivers that are not in SSM and have weird biogeochemistry values
# (doesn't include overlapping rivers between SSM and LO)
weird_rivers = ['Brooks Peninsula', 'Campbell River', 'Clayoquot', 'Holberg',
                'Homathco River', 'Klinaklini River', 'Knight Inlet',
                'Neil Creek', 'Nimpkish River', 'North East Vancouver Island',
                'Owikeno Lake', 'Salmon River', 'Seymour Inlet', 'Tahsis',
                'Toba Inlet', 'Tsitika River', 'Willamette R']

# loop through all rivers
for i,rname in enumerate(rivnames):

    print('{}: {}'.format(i,rname))

    # get river index
    rID = rivids[i]
    
    riv_fn = ''

    # find Ecology's timeseries based on river id
    for fn in all_his_fns:
        root, ext = os.path.splitext(fn)
        if root.startswith(str(rID)) and ext == '.xlsx':
            riv_fn = fn
    
    # Let user know if couldn't find timeseries for a given river
    if riv_fn == '0':
        print('No history file found for {}'.format(rname))
        continue
    # Otherwise, read history file as df
    else:
        riv_fp = str(riv_dir)  + '/' + riv_fn
        riv_df = pd.read_excel(riv_fp, skiprows=[0])

    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    riv_df = riv_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)

    # replace all zeros with nans, so zeros don't bias data
    riv_df = riv_df.replace(0, np.nan)

    # calculate averages
    # (compress 1999-2017 timeseries to single day, with an average for each day)
    riv_avgs_df = riv_df.groupby(['Month','Day']).mean().reset_index()
    # calculate standard deviation
    riv_sds_df = riv_df.groupby(['Month','Day']).std().reset_index()

    # replace all nans with zeros, so I'm no longer injecting nans
    riv_avgs_df = riv_avgs_df.replace(np.nan,0)
    riv_sds_df = riv_sds_df.replace(np.nan,0)

    # Plot and save averages for each source
    fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
    ax = axes.ravel()
    for j,vn in enumerate(vns):
        i = j+1
        # label subplot
        ax[i].set_title(vn,fontsize=14)
        # Plot individual years
        for yr in range(1999,2017):
            riv_yr_df = riv_df.loc[riv_df['Year'] == yr]
            # Insert a nan on Feb 29 if not a leap year
            if np.mod(yr,4) != 0:
                nans = [np.nan]*29
                riv_yr_df = riv_yr_df.reset_index(drop=True) # reset all dataframes to index from 0
                riv_yr_df.loc[58.5] = nans # leap year is 60th Julian day, so add a new 59th index since Python indexes from 0
                riv_yr_df = riv_yr_df.sort_index().reset_index(drop=True) # sort indices and renumber
            if yr == 2017:
                yrday_17 = pd.date_range(start ='1/1/2020', end ='8/02/2020', freq ='D') # don't have full 2017 dataset
                ax[i].plot(yrday_17,riv_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
            else:
                ax[i].plot(yrday,riv_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
        # Plot average
        ax[i].plot(yrday,riv_avgs_df[vn].values, label='average', color='black', linewidth=1.5)
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
    plt.suptitle(rname,fontsize=18)
    # Save figure if not a weird river
    if rname not in weird_rivers:
        figname = rname + '.png'
        save_path = clim_dir / 'climatology_plots' / figname
        fig.savefig(save_path)
        plt.close('all')
        # plt.show()
    
    # Add data to climatology dataframes, and convert to units that LiveOcean expects
    flow_clim_df = pd.concat([flow_clim_df, pd.Series(riv_avgs_df['Flow(m3/s)'].values, name=rname)], axis = 1)     # [m3/s]
    temp_clim_df = pd.concat([temp_clim_df, pd.Series(riv_avgs_df['Temp(C)'].values, name=rname)], axis = 1)        # [C]
    # don't add biogeochem for weird rivers with unrealistic values (pad with nans)
    if rname in weird_rivers:
        nans = np.empty((366,))
        nans[:] = np.nan
        NO3_clim_df  = pd.concat([NO3_clim_df, pd.Series(nans, name=rname)], axis = 1)
        NH4_clim_df  = pd.concat([NH4_clim_df, pd.Series(nans, name=rname)], axis = 1)
        TIC_clim_df  = pd.concat([TIC_clim_df, pd.Series(nans, name=rname)], axis = 1)
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(nans, name=rname)], axis = 1)
        DO_clim_df   = pd.concat([DO_clim_df, pd.Series(nans, name=rname)], axis = 1)
    # add biogeochem for all normal rivers
    else:
        NO3_clim_df  = pd.concat([NO3_clim_df, pd.Series(riv_avgs_df['NO3+NO2(mg/L)'] * 71.4, name=rname)], axis = 1)  # [mmol/m3]
        NH4_clim_df  = pd.concat([NH4_clim_df, pd.Series(riv_avgs_df['NH4(mg/L)'] * 71.4, name=rname)], axis = 1)      # [mmol/m3]
        TIC_clim_df  = pd.concat([TIC_clim_df, pd.Series(riv_avgs_df['DIC(mmol/m3)'], name=rname)], axis = 1)          # [mmol/m3]
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(riv_avgs_df['Alk(mmol/m3)'], name=rname)], axis = 1)          # [meq/m3]
        DO_clim_df   = pd.concat([DO_clim_df, pd.Series(riv_avgs_df['DO(mg/L)'] * 31.26, name=rname)], axis = 1)      # [mmol/m3]

# Calculate summary statistics for all tinyrivs
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

# fill weird river biogeochemisty values with average values
for rname in weird_rivers:
    # temp_clim_df[rname] = clim_avgs['Temp(C)']                # [C]
    NO3_clim_df[rname]  = clim_avgs['NO3+NO2(mg/L)'] * 71.4   # [mmol/m3]
    NH4_clim_df[rname]  = clim_avgs['NH4(mg/L)'] * 71.4       # [mmol/m3]
    TIC_clim_df[rname]  = clim_avgs['DIC(mmol/m3)']           # [mmol/m3]
    Talk_clim_df[rname] = clim_avgs['Alk(mmol/m3)']           # [meq/m3]
    DO_clim_df[rname]   = clim_avgs['DO(mg/L)'] * 31.26       # [mmol/m3]
    # Plot and save averages for weird rivers
    fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
    ax = axes.ravel()
    for j,vn in enumerate(vns):
        i = j+1
        # label subplot
        ax[i].set_title(vn,fontsize=14)
        # Plot actual river flowrate
        if vn == 'Flow(m3/s)' or vn == 'Temp(C)':
            ax[i].plot(yrday,clim_df_list[j][rname].values, label='Climatology profile for this river', color='black', linewidth=1.5)
        else:
            scale = 1
            if vn == 'NO3+NO2(mg/L)' or vn == 'NH4(mg/L)':
                scale = 71.4
            elif vn == 'DO(mg/L)':
                scale = 31.26
            # Plot average river climatology for biogeochemistry
            ax[i].plot(yrday,clim_df_list[j][rname].values/scale, label='Climatology profile for this river', color='black', linewidth=1.5) # this line here just for legend
            ax[i].plot(yrday,clim_df_list[j][rname].values/scale, label='Filled with average climatology of all rivers', color='silver', linewidth=1.5)
        # fontsize of tick labels
        ax[i].tick_params(axis='both', which='major', labelsize=12)
        ax[i].tick_params(axis='x', which='major', rotation=30)
        ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
        # create legend
        if i ==7:
            handles, labels = ax[7].get_legend_handles_labels()
            ax[0].legend(handles, labels, loc='center', ncol = 1,fontsize=14)
            ax[0].axis('off')
        # Define the date format
        if i >= 6:
            date_form = mdates.DateFormatter("%b")
            ax[i].xaxis.set_major_formatter(date_form)
    # plot title is name of source
    plt.suptitle(rname,fontsize=18)
    figname = rname + '.png'
    save_path = clim_dir / 'climatology_plots' / figname
    fig.savefig(save_path)
    plt.close('all')
    # plt.show()

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
plt.suptitle('River Climatology Summary (n={})'.format(len(rivnames)),fontsize=18)
# Save figure
figname = 'river_summary.png'
save_path = clim_dir / figname
fig.savefig(save_path)
# plt.close('all')
plt.show()