"""
Make climatologies for Mohamedali et al. (2020) tiny rivers.
Discharge rate, temperature, and biogeochemisty variables.

Based on Mohamedali et al. (2020) timeseries, using data stored in 
LO_data/trapsD01/processed_data/river_data_mohamedali_etal_2020.nc

        To run, from ipython:
        run make_climatology_moh20_tinyrivs.py

To create individual climatology figures, run from ipython with:
run make_climatology_moh20_tinyrivs.py -test True

Figures saved in:
LO_output/pre/trapsP01/tiny_rivers/[ctag]/Data_historical/climatology_plots

Note that running with -test True adds
a couple minutes to run time. (~ 2 min)
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
clim_dir = Ldir['LOo'] / 'pre' / trapsP_vers / 'moh20_tinyrivers' / ctag / 'Data_historical'
Lfun.make_dir(clim_dir)

# get flow and loading data
moh20_triv_fn = Ldir['data'] / trapsD / 'processed_data'/ 'river_data_mohamedali_etal_2020.nc'
moh20_river_data_ds = xr.open_dataset(moh20_triv_fn)

# get triv names
trivnames_all = moh20_river_data_ds['name'].values

# Remove pre-existing LO rivers
# read overlapping rivers
repeatrivs_fn = Ldir['data'] / trapsD / 'LiveOcean_SSM_rivers.xlsx'
repeatrivs_df = pd.read_excel(repeatrivs_fn)
SSM_repeats = repeatrivs_df['SSM_rname'].tolist()
# remove nans
SSM_repeats = [x for x in SSM_repeats if str(x) != 'nan']
# remove repeat river names from list of river names & remove Willamette River
SSM_repeats_and_Willamette = SSM_repeats + ['Willamette R']
trivnames = [river for river in trivnames_all if river not in SSM_repeats_and_Willamette]

# # just test 5 trivs for now -------------------------------------------------
# trivnames = trivnames[28:33]

# separate rivers with weird biogeochem
weird_biogeochem = ['Neil Creek', 'Seymour Inlet', 'Holberg',
                    'Owikeno Lake', 'Salmon River', 'Brooks Peninsula',
                    'Clayoquot', 'Toba Inlet', 'Homathco River ',
                    'Campbell River', 'Tsitika River', 'Nimpkish River',
                    'Tahsis', 'Alberni Inlet', 'Knight Inlet',
                    'Willamette R', 'Klinaklini River',
                    'North East Vancouver Island']
weird_DO = ['Victoria_SJdF']
weird_DO_and_NH4 = ['Vancouver Isl C']
# get list of rivers with more realistic biogeochem
trivnames_realisticBGC = [river for river in trivnames
                          if  river not in weird_biogeochem
                          and river not in weird_DO
                          and river not in weird_DO_and_NH4]

# initialize dataframes for all trivs
DO_clim_df   = pd.DataFrame()
flow_clim_df = pd.DataFrame()
temp_clim_df = pd.DataFrame()
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
TIC_clim_df  = pd.DataFrame()
Talk_clim_df = pd.DataFrame()

# variable names
# IMPORTANT: all variables must be in the same order in the following lests
# TO-DO: turn this into a dictionay so order will not matter
print('TO-DO: MAKE A DICTIONARY OF VARIABLE NAMES (including clims later in code)')
vns = ['Flow(m3/s)','NO3(mmol/m3)','NH4(mmol/m3)',
       'TIC(mmol/m3)','Talk(meq/m3)','DO(mmol/m3)','Temp(C)']
pickle_names = ['flow', 'NO3', 'NH4',
                 'TIC',  'Talk', 'DO', 'temp']
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']

# create one-year date range for plotting
yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')

# Replace 2012-2014 data with nans (so we don't bias climatologies)
# Do this for rivers in which those data were accidentally shifted by 3 months
shifted = ['Kitsap NE', 'Kitsap_Hood', 'Lynch Cove',
           'NW Hood', 'Port Gamble', 'Tahuya']
keys = ['flow', 'temp', 'NO3', 'NH4', 'TIC', 'Talk', 'DO']
# create a mask for the years 2012-2014
mask_2012_2014 = (moh20_river_data_ds['date'] > np.datetime64('2012-07-31')) & (moh20_river_data_ds['date.year'] < 2015)

# replace the 2012-2014 flow data with nan for the shifted years
for riv in shifted:
    # get river index
    riv_index = int(moh20_river_data_ds.where(moh20_river_data_ds['name'] == riv, drop=True)['source'])
    # replace 2012-2014 flow data with nans
    for key in keys:
        moh20_river_data_ds[key].loc[dict(source=riv_index, date=mask_2012_2014)] = np.nan

#################################################################################
#                          Calculate climatologies                              #
#################################################################################

print('Calculating tiny river climatologies...')

# loop through all nonpoint sources
for i,rname in enumerate(trivnames):

    # turn dataset information for this wwtp into a dataframe
    triv_df, triv_avgs_df, triv_sds_df = traps_helper.ds_to_avgdf(rname,moh20_river_data_ds)

#################################################################################
#                            Create climatologies                               #
#################################################################################

    # Add data to climatology dataframes
    flow_clim_df = pd.concat([flow_clim_df, pd.Series(triv_avgs_df['Flow(m3/s)'].values, name=rname)], axis = 1)# [m3/s]
    temp_clim_df = pd.concat([temp_clim_df, pd.Series(triv_avgs_df['Temp(C)'].values, name=rname)], axis = 1)   # [C]
    # create list of nans
    nans = np.empty((366,))
    nans[:] = np.nan
    # don't add biogeochem for weird rivers with unrealistic values (pad with nans)
    if rname in weird_biogeochem:
        NO3_clim_df  = pd.concat([NO3_clim_df,  pd.Series(nans, name=rname)], axis = 1)
        NH4_clim_df  = pd.concat([NH4_clim_df,  pd.Series(nans, name=rname)], axis = 1)
        TIC_clim_df  = pd.concat([TIC_clim_df,  pd.Series(nans, name=rname)], axis = 1)
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(nans, name=rname)], axis = 1)
        DO_clim_df   = pd.concat([DO_clim_df,   pd.Series(nans, name=rname)], axis = 1)
    elif rname in weird_DO:
        DO_clim_df   = pd.concat([DO_clim_df,   pd.Series(nans, name=rname)], axis = 1)
        NO3_clim_df  = pd.concat([NO3_clim_df,  pd.Series(triv_avgs_df['NO3(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        NH4_clim_df  = pd.concat([NH4_clim_df,  pd.Series(triv_avgs_df['NH4(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        TIC_clim_df  = pd.concat([TIC_clim_df,  pd.Series(triv_avgs_df['TIC(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(triv_avgs_df['Talk(meq/m3)'], name=rname)], axis = 1) # [meq/m3]
    elif rname in weird_DO_and_NH4:
        DO_clim_df   = pd.concat([DO_clim_df,   pd.Series(nans, name=rname)], axis = 1)
        NH4_clim_df  = pd.concat([NH4_clim_df,  pd.Series(nans, name=rname)], axis = 1)
        NO3_clim_df  = pd.concat([NO3_clim_df,  pd.Series(triv_avgs_df['NO3(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        TIC_clim_df  = pd.concat([TIC_clim_df,  pd.Series(triv_avgs_df['TIC(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(triv_avgs_df['Talk(meq/m3)'], name=rname)], axis = 1) # [meq/m3]
    # add biogeochem for all normal rivers
    else:
        NO3_clim_df  = pd.concat([NO3_clim_df,  pd.Series(triv_avgs_df['NO3(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        NH4_clim_df  = pd.concat([NH4_clim_df,  pd.Series(triv_avgs_df['NH4(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        TIC_clim_df  = pd.concat([TIC_clim_df,  pd.Series(triv_avgs_df['TIC(mmol/m3)'], name=rname)], axis = 1) # [mmol/m3]
        Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(triv_avgs_df['Talk(meq/m3)'], name=rname)], axis = 1) # [meq/m3]
        DO_clim_df   = pd.concat([DO_clim_df,   pd.Series(triv_avgs_df['DO(mmol/m3)'],  name=rname)], axis = 1) # [mmol/m3]

# Calculate summary statistics for all trivs
clim_avgs = pd.DataFrame()
clim_max = pd.DataFrame()
clim_min = pd.DataFrame()
clim_sds = pd.DataFrame()

# climatology dfs
clims = [flow_clim_df, NO3_clim_df, NH4_clim_df,
        TIC_clim_df, Talk_clim_df, DO_clim_df, temp_clim_df]

for i,vn in enumerate(vns):
    # average values of all nonpoint sources
    clim_avgs[vn] = clims[i].mean(axis=1)
    # max climatology values
    clim_max[vn] = clims[i].max(axis=1)
    # min climatology values
    clim_min[vn] = clims[i].min(axis=1)
    # standard deviation of all nonpoint sources
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
plt.suptitle('Tiny River Climatology Summary (n={})'.format(len(trivnames_realisticBGC)),fontsize=18)
# Save figure
figname = 'tiny_river_summary.png'
save_path = clim_dir / figname
fig.savefig(save_path)

#################################################################################
# Create climatology for rivers w/ weird biogeochem (using avg of other rivers) #
#################################################################################

# TODO: add average climatology based on watershed,
# so we don't bias rivers in urban and remote environments

# fill weird river biogeochemisty values with average values
for rname in weird_biogeochem:
    NO3_clim_df[rname]  = clim_avgs['NO3(mmol/m3)']     # [mmol/m3]
    NH4_clim_df[rname]  = clim_avgs['NH4(mmol/m3)']     # [mmol/m3]
    TIC_clim_df[rname]  = clim_avgs['TIC(mmol/m3)']     # [mmol/m3]
    Talk_clim_df[rname] = clim_avgs['Talk(meq/m3)']     # [meq/m3]
    DO_clim_df[rname]   = clim_avgs['DO(mmol/m3)']      # [mmol/m3]

# fill weird river DO values with average values
for rname in weird_DO:
    DO_clim_df[rname]   = clim_avgs['DO(mmol/m3)']      # [mmol/m3]

# fill weird river DO + NH4 values with average values
for rname in weird_DO_and_NH4:
    NH4_clim_df[rname]  = clim_avgs['NH4(mmol/m3)']     # [mmol/m3]
    DO_clim_df[rname]   = clim_avgs['DO(mmol/m3)']      # [mmol/m3]

print('Tiny river climatologies done\n')

#################################################################################
#                Plot all individual tiny river climatologies                   #
#################################################################################

if args.testing == True:
    print('Plotting...')

    # generate directory to save files
    fig_dir = clim_dir / 'climatology_plots'
    Lfun.make_dir(fig_dir)

    for j,rname in enumerate(trivnames):

        print('{}/{}: {}'.format(j+1,len(trivnames),rname))

        # turn dataset information for this triv into a dataframe
        triv_df, triv_avgs_df, triv_sds_df = traps_helper.ds_to_avgdf(rname,moh20_river_data_ds)

        # Plot climatologies for each source
        fig, axes = plt.subplots(7,1, figsize=(10, 9), sharex=True)
        ax = axes.ravel()

        # loop through variables
        for i,vn in enumerate(vns):

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
            values_to_plot = triv_df[vn].values*scale # scale variable properly
            time = triv_df['Date'].values
            ax[i].plot(time,values_to_plot, color='deeppink', label='Raw data',
                       linewidth=2, alpha=0.5)

            # Plot average
            leapyear_clim = clims[i][rname].values*scale #triv_avgs_df[vn].values*scale
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
            if i > 0 and (rname in weird_biogeochem or rname in weird_DO or rname in weird_DO_and_NH4):
                maxval = np.nanmax(leapyear_clim)
            else:
                maxval = 1.3*np.nanmax(values_to_plot)
            ax[i].set_ylim([0,maxval])

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
        plt.suptitle(rname,fontsize=14)
        plt.subplots_adjust(top=0.95)

        # Save figure
        figname = rname + '.png'
        save_path = clim_dir / 'climatology_plots' / figname
        fig.savefig(save_path)
        plt.close('all')

    print('Done')

#################################################################################
#                             Save climatologies                                #
#################################################################################

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