"""
Script to plot and compare new WWTP data and my older climatologies.
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cmocean
import datetime as dt
import matplotlib.dates as mdates
from pathlib import Path

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

#################################################################################
#                             Facility information                              #
#################################################################################

plt.close('all')
# where to put output figures
out_dir = Ldir['LOo'] / 'loading_test' / 'point_source_integration'
Lfun.make_dir(out_dir / 'Wasielewski_etal2024' / 'WWTPs')
Lfun.make_dir(out_dir / 'Mohamedali_etal2020' / 'WWTPs')

print('\n')

#################################################################################
#                   Wasielewski et al., 2024 loading plots                      #
#################################################################################

plot_individual_WWTP_loads = False
plot_facility_type_load_comparison = True

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent.parent.parent
with open(this_dir /'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# location of point source data to process
pointsource_dir = Ldir['data'] / trapsD / 'wasielewski_etal2024' / 'point_sources'
pointsource_meta = pointsource_dir / 'fac_attributes.csv'
pointsource_loads = pointsource_dir / 'nutrient_loads.csv'
# get point source data
psmeta_df = pd.read_csv(pointsource_meta)
psloads_df = pd.read_csv(pointsource_loads)

facilities_all = psmeta_df['FAC_ID'].values.tolist()

# separate by facility type
WWTP_df = psmeta_df[psmeta_df['FAC_TYPE']=='sic_4952']
hatcheries_df = psmeta_df[psmeta_df['FAC_TYPE']=='sic_0921']
industrial_df = psmeta_df[psmeta_df['FAC_TYPE']=='sic_INDU']
              
# Plot individual loading plots for WWTPs -----------------------------------------------------
if plot_individual_WWTP_loads == True:
        print('Plotting all WWTP loads (Wasielewski et al., 2024)...')
        for fac_ID in WWTP_df['FAC_ID']:

                # get WWTP ID
                fac_ID = fac_ID

                # get name
                fac_name = WWTP_df.loc[WWTP_df['FAC_ID'] == fac_ID, 'FAC_NAME'].values[0]

                # get flow corresponding to desired WWTP
                # replace any '.' with NaN
                no_empties_flow = psloads_df.loc[psloads_df['FAC_ID'] == fac_ID, 'FLOW_MGD'].replace('.', np.nan).astype(float)
                new_flow_millGall_day = no_empties_flow.values
                # convert to m3/s
                new_flow_m3_s = new_flow_millGall_day * 0.0438126364

                # get nutrient data [mg/L]
                no_empties_NO3 = psloads_df.loc[psloads_df['FAC_ID'] == fac_ID, 'NO2NO3N_MG_L'].replace('.', np.nan).astype(float)
                no_empties_NH4 = psloads_df.loc[psloads_df['FAC_ID'] == fac_ID, 'NH4N_MG_L'].replace('.', np.nan).astype(float)
                new_NO3 = no_empties_NO3.values
                new_NH4 = no_empties_NH4.values
                # sum NO3 and NH4
                new_nutrients = new_NO3 + new_NH4 # [mg/L]

                # get total nutrient loading [kg/day]
                new_nutrient_load_m3mg_sL = new_flow_m3_s * new_nutrients # [m3/s * mg/L]
                # convert to kg/day: [1 kg = 1/1000 m3mg/L] [1day = 60*60*24sec]
                new_nutrient_load_kg_day = new_nutrient_load_m3mg_sL * 60 * 60 * 24 / 1000

                # rename for consistency
                new_nutrient_load_15years = new_nutrient_load_kg_day

                ## Plot ----------------

                # get time array
                month = psloads_df.loc[psloads_df['FAC_ID'] == fac_ID, 'MONTH'].replace('.', np.nan).astype(int)
                year = psloads_df.loc[psloads_df['FAC_ID'] == fac_ID, 'YEAR'].replace('.', np.nan).astype(int)
                # create time array
                t_new = np.array([dt.datetime(year, month, 1) for year, month in zip(year, month)], dtype='datetime64[M]') 

                plt.close('all')
                fig, ax = plt.subplots(1,1,figsize = (13,7))

                # plot new data
                ax.plot(t_new,new_nutrient_load_15years,marker='o',
                        linewidth=1,color='hotpink', alpha=0.8)

                # format figure
                ax.set_title('{} nutrient load'.format(fac_name),
                        fontsize=14,fontweight='bold')
                ax.set_ylabel('Nutrient load [kg/day]',fontsize=12)

                ax.set_ylim([0,np.nanmax(new_nutrient_load_15years)*1.1])
                ax.set_xlim([pd.to_datetime('2005-01-01'),pd.to_datetime('2020-12-31')])
                # ax.set_xlim([pd.to_datetime('2013-01-01'),pd.to_datetime('2017-12-31')])
                ax.xaxis.set_major_locator(mdates.YearLocator())
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
                ax.tick_params(axis='x', labelrotation=30)
                ax.tick_params(axis='both', labelsize=12)
                ax.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')

                plt.savefig(out_dir / 'Wasielewski_etal2024' / 'WWTPs' / (fac_ID +'.png'))
        print('    Done')

# Create nutrient load summary of the different types of plants -------------------------------
if plot_facility_type_load_comparison == True:
        
        fig, ax = plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]},figsize = (7,9))

        # plot outfall locations on a map
        # get the grid data
        ds = xr.open_dataset('../../../../../LO_data/grids/cas7/grid.nc')
        z = -ds.h.values
        mask_rho = np.transpose(ds.mask_rho.values)
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        X = lon[0,:] # grid cell X values
        Y = lat[:,0] # grid cell Y values
        plon, plat = pfun.get_plon_plat(lon,lat)
        # make a version of z with nans where masked
        zm = z.copy()
        zm[np.transpose(mask_rho) == 0] = np.nan
        zm[np.transpose(mask_rho) != 0] = -1
        # add land and water mask to both subplots
        ax[0].pcolormesh(plon, plat, zm, vmin=-20, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
        # add WWTPs
        ax[0].scatter(WWTP_df['FAC_LON'],WWTP_df['FAC_LAT'], s=30,zorder=5,
                facecolors='deeppink', edgecolors='none', alpha=0.9,
                label='WWTPs')
        # add hatcheries
        ax[0].scatter(hatcheries_df['FAC_LON'],hatcheries_df['FAC_LAT'], s=30, zorder=5,
                facecolors='black', edgecolors='none', alpha=0.9, marker='d',
                label='Hatcheries')
        # add industrial facilities
        ax[0].scatter(industrial_df['FAC_LON'],industrial_df['FAC_LAT'], s=30, zorder=5,
                facecolors='deepskyblue', edgecolors='none', alpha=0.9, marker='s',
                label='Industrial Facilities')
        # format figure
        pfun.dar(ax[0])
        pfun.add_coast(ax[0],color='paleturquoise')
        ax[0].set_title('Facility locations',fontsize=14)
        ax[0].set_ylim([46.8,49.6])
        ax[0].set_xlim([-125,-121])
        ax[0].legend(loc='upper left',fontsize=12)

        # plot nutrient loads (monthly)
        # get all loading data for each facility type
        WWTP_load_df = psloads_df[psloads_df['FAC_ID'].isin(WWTP_df['FAC_ID'])]
        hatchery_load_df = psloads_df[psloads_df['FAC_ID'].isin(hatcheries_df['FAC_ID'])]
        industrial_load_df = psloads_df[psloads_df['FAC_ID'].isin(industrial_df['FAC_ID'])]
        # first, sum all of the individual plants together to get one 2005 - 2020 timeseries (which is the sum of all plants)
        all_WWWTP_loads = WWTP_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        all_hatchery_loads = hatchery_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        all_industrial_loads = industrial_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        # then, get the mean annual loading profile (still, sum of all plants)
        annualmean_all_WWTP_loads = all_WWWTP_loads.groupby(['MONTH']).mean()
        annualmean_all_hatchery_loads = all_hatchery_loads.groupby(['MONTH']).mean()
        annualmean_all_industrial_loads = all_industrial_loads.groupby(['MONTH']).mean()
        # plot time series of loads
        t = pd.date_range(start='2005-01-01', end='2005-12-01', freq='MS')
        month_days = [1/31, 1/28.25, 1/31, 1/30,
                      1/31, 1/30,    1/31, 1/31,
                      1/30, 1/31,    1/30, 1/31]
        ax[1].plot(t,annualmean_all_WWTP_loads*month_days, color='deeppink',linewidth=3,alpha=0.7)
        ax[1].plot(t,annualmean_all_hatchery_loads*month_days, color='black',linewidth=3,alpha=0.7)
        ax[1].plot(t,annualmean_all_industrial_loads*month_days, color='deepskyblue',linewidth=3,alpha=0.7)
        # format figure
        ax[1].set_title('Mean Annual Nutrient Loads',fontsize=14)
        ax[1].set_ylabel('Total Nitrogen Load\n[kg/day]',fontsize=12)
        ax[1].set_yscale('log')
        ax[1].set_ylim([1,1e6])
        ax[1].set_xlim([t[0],t[-1]])
        ax[1].xaxis.set_major_locator(mdates.MonthLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

        # format and save figure
        plt.suptitle('Wasielewski et al. (2024)',fontsize=14,fontweight='bold')
        plt.tight_layout
        plt.show()
        plt.savefig(out_dir / 'Wasielewski_etal2024' / 'all_loads_comparison.png')

# Save WWTP data for later analysis
wasielewski24_WWTP_df = psloads_df.loc[psloads_df['FAC_ID'].isin(WWTP_df['FAC_ID'])]

#################################################################################
#                    Mohamedali et al., 2020 loading plots                      #
#################################################################################

plot_individual_WWTP_loads_mohamedali = False
plot_facility_type_load_comparison_mohamedali = True

# location of point source data to process
pointsource_dir = Ldir['data'] / trapsD / 'mohamedali_etal2020'
pointsource_loads = pointsource_dir / 'all_point_source_data.nc'
# get point source data
psloads_ds = xr.open_dataset(pointsource_loads)

# get WWTPs vs. industrial facilities
industrial_names = ['BP Cherry Point',
                    'Conoco Phillips',
                    'Intalco',
                    'Kimberly_Clark',
                    'Nippon Paper',
                    'Port Townsend Paper',
                    'Shell Oil',
                    'Tesoro',
                    'US Oil & Refining',
                    'West Rock']
WWTP_names = [fac for fac in psloads_ds['name'].values if fac not in industrial_names]

# Plot individual loading plots for WWTPs -----------------------------------------------------
if plot_individual_WWTP_loads_mohamedali == True:
        print('Plotting all WWTP loads (Mohamedali et al., 2020)...')

        for fac_name in WWTP_names:

                # get flow [m3/s] & subsample to one value per month
                raw_flow = psloads_ds.sel(source=(psloads_ds['name']==fac_name))['flow'].resample(date='1MS').first()[0]

                # get NO3 and NH4 [mmol/m3] & subsample to one value per month
                raw_NO3 = psloads_ds.sel(source=(psloads_ds['name']==fac_name))['NO3'].resample(date='1MS').first()[0]
                raw_NH4 = psloads_ds.sel(source=(psloads_ds['name']==fac_name))['NH4'].resample(date='1MS').first()[0]
                # sum nutrients
                raw_nutrients = raw_NO3 + raw_NH4 # [mmol/m3]

                # get total nutrient loading [kg/day]
                raw_nutrient_load_mmol_s = raw_flow * raw_nutrients # [m3/s * mmol/m3 = mmol/s]
                # convert to kg/day: [1g = 14.01/1000^2mmol] [1day = 60*60*24sec]
                raw_nutrient_load_kg_day = raw_nutrient_load_mmol_s * 14.01 * 60 * 60 * 24 / 1000 / 1000

                # Plot ----------------

                # create time array
                t = psloads_ds.date.resample(date='1MS').first()
                # initialize figure
                fig, ax = plt.subplots(1,1,figsize = (13,7))
                # plot loads
                ax.plot(t,raw_nutrient_load_kg_day,marker='o',
                        linewidth=1,color='hotpink', alpha=0.8)
                # format figure
                ax.set_title('{} nutrient load'.format(fac_name),
                        fontsize=14,fontweight='bold')
                ax.set_ylabel('Nutrient load [kg/day]',fontsize=12)

                ax.set_ylim([0,np.nanmax(raw_nutrient_load_kg_day)*1.1])
                ax.set_xlim([pd.to_datetime('1999-01-01'),pd.to_datetime('2017-12-31')])
                ax.xaxis.set_major_locator(mdates.YearLocator())
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
                ax.tick_params(axis='x', labelrotation=30)
                ax.tick_params(axis='both', labelsize=12)
                ax.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')

                plt.savefig(out_dir / 'Mohamedali_etal2020' / 'WWTPs' / (fac_name +'.png'))
                plt.close()
        print('    Done')

# Create nutrient load summary of the different types of plants -------------------------------
if plot_facility_type_load_comparison_mohamedali == True:
        
        fig, ax = plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]},figsize = (7,9))

        # plot outfall locations on a map
        # add land and water mask to both subplots
        ax[0].pcolormesh(plon, plat, zm, vmin=-20, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
        # add WWTPs
        ax[0].scatter(psloads_ds.sel(source=psloads_ds['name'].isin(WWTP_names))['lon'],
                psloads_ds.sel(source=psloads_ds['name'].isin(WWTP_names))['lat'],
                s=30,zorder=5,facecolors='deeppink', edgecolors='none', alpha=0.9,
                label='WWTPs')
        # add industrial facilities
        ax[0].scatter(psloads_ds.sel(source=psloads_ds['name'].isin(industrial_names))['lon'],
                psloads_ds.sel(source=psloads_ds['name'].isin(industrial_names))['lat'],
                s=30, zorder=5,facecolors='deepskyblue', edgecolors='none', alpha=0.9, marker='s',
                label='Industrial Facilities')
        # format figure
        pfun.dar(ax[0])
        pfun.add_coast(ax[0],color='paleturquoise')
        ax[0].set_title('Facility locations',fontsize=14)
        ax[0].set_ylim([46.8,49.6])
        ax[0].set_xlim([-125,-121])
        ax[0].legend(loc='lower left',fontsize=12)

        # plot nutrient loads (monthly)
        # create new variable representing nutrient load in kg/day
        psloads_ds['load [kg/d]'] = psloads_ds['flow'] * (psloads_ds['NO3'] + psloads_ds['NH4']) * 14.01 * 60 * 60 * 24 / 1000 / 1000
        # get all loading data for each facility type
        WWTP_load_ds = psloads_ds.sel(source=psloads_ds['name'].isin(WWTP_names))
        industrial_load_ds = psloads_ds.sel(source=psloads_ds['name'].isin(industrial_names))
        # get sum of load from all sources
        all_WWTP_loads = WWTP_load_ds['load [kg/d]'].sum(dim='source')
        all_industrial_loads = industrial_load_ds['load [kg/d]'].sum(dim='source')
        # then, get the mean annual loading profile with monthly resolution (still, sum of all plants)
        annualmean_all_WWTP_loads = all_WWTP_loads.groupby('date.month').mean(dim='date')
        annualmean_all_industrial_loads = all_industrial_loads.groupby('date.month').mean(dim='date')
        # plot time series of loads
        t = pd.date_range(start='2005-01-01', end='2005-12-01', freq='MS')
        ax[1].plot(t,annualmean_all_WWTP_loads, color='deeppink',linewidth=3,alpha=0.7)
        ax[1].plot(t,annualmean_all_industrial_loads, color='deepskyblue',linewidth=3,alpha=0.7)

        # format figure
        ax[1].set_title('Mean Annual Nutrient Loads',fontsize=14)
        ax[1].set_ylabel('Total Nitrogen Load\n[kg/day]',fontsize=12)
        ax[1].set_yscale('log')
        ax[1].set_ylim([1,1e6])
        ax[1].set_xlim([t[0],t[-1]])
        ax[1].xaxis.set_major_locator(mdates.MonthLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

        # format and save figure
        plt.suptitle('Mohamedali et al. (2020)',fontsize=14,fontweight='bold')
        plt.tight_layout
        plt.show()
        plt.savefig(out_dir / 'Mohamedali_etal2020' / 'all_loads_comparison.png')

# Save WWTP data for later analysis
mohamedali2020_WWTP_ds = psloads_ds.sel(source=psloads_ds['name'].isin(WWTP_names))

#################################################################################
#                             Dataset comparisons                               #
#################################################################################

compare_plant_locations = True
largest_plants_analysis = True
unique_plant_locations = True

# plot WWTP outfall locations on an interactive map ----------------------------------
if compare_plant_locations == True:
        fig, ax = plt.subplots(1,1,figsize = (12,9))
        # get the grid data
        ds = xr.open_dataset('../../../../../LO_data/grids/cas7/grid.nc')
        z = -ds.h.values
        mask_rho = np.transpose(ds.mask_rho.values)
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        X = lon[0,:] # grid cell X values
        Y = lat[:,0] # grid cell Y values
        plon, plat = pfun.get_plon_plat(lon,lat)
        # make a version of z with nans where masked
        zm = z.copy()
        zm[np.transpose(mask_rho) == 0] = np.nan
        zm[np.transpose(mask_rho) != 0] = -1
        # add land and water mask to both subplots
        ax.pcolormesh(plon, plat, zm, vmin=-20, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # add old WWTP data
        lats = mohamedali2020_WWTP_ds['lat'].values
        lons = mohamedali2020_WWTP_ds['lon'].values
        ax.scatter(lons,lats, color='deeppink',s=30,alpha=0.8,
                zorder=5,label='Mohamedali et al., 2020')

        # add new WWTP data
        lats = WWTP_df['FAC_LAT'].values
        lons = WWTP_df['FAC_LON'].values
        ax.scatter(lons,lats, s=30,zorder=5,facecolors='none',
                edgecolors='black',
                label='Wasielewski et al., 2024')
        
        # show labels
        for i,wwtp in enumerate(mohamedali2020_WWTP_ds['name'].values):
                ax.text(mohamedali2020_WWTP_ds['lon'].values[i],
                        mohamedali2020_WWTP_ds['lat'].values[i],
                        wwtp, color='deeppink',
                        horizontalalignment='right')
        for i,wwtp in enumerate(WWTP_df['FAC_NAME'].values):
                ax.text(WWTP_df['FAC_LON'].values[i],
                        WWTP_df['FAC_LAT'].values[i],
                        wwtp, color='k')

        # format figure
        pfun.dar(ax)
        pfun.add_coast(ax,color='paleturquoise')
        ax.set_title('WWTP locations\npan and zoom to locate overlapping WWTPs',fontsize=14)
        ax.set_ylim([46.8,49.4])
        ax.set_xlim([-125,-121])
        ax.legend(loc='lower left',fontsize=14)

        plt.show()

# Get the largest plants in each dataset ---------------------------------
if largest_plants_analysis == True: 
       
        # create dictionaries of WWTP:avg_load[kg/day] for both datasets
        # initialize dicts
        wasielewski24_load_dict = {}
        mohamedali2020_load_dict = {}

        # add values to dictionary

        # Wasielewski et al., 2024
        for ID in wasielewski24_WWTP_df['FAC_ID'].unique():
                # get wwtp name
                wwtp = psmeta_df.loc[psmeta_df['FAC_ID'] == ID, 'FAC_NAME'].values[0]
                # get mean loading of each WWTP
                flow = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'] == ID, 'FLOW_MGD'].replace('.', np.nan).astype(float) * 0.0438126364 # m3/s from mill gallons per day
                NO3 = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'] == ID, 'NO2NO3N_MG_L'].replace('.', np.nan).astype(float) # mg/L
                NH4 = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'] == ID, 'NH4N_MG_L'].replace('.', np.nan).astype(float) # mg/L
                # calculate loads
                load = (flow * (NO3+NH4) * 60 * 60 * 24 / 1000).values # convert to kg/day: [1 kg = 1/1000 m3mg/L] [1day = 60*60*24sec]
                # replace all zeros with nans, so that they don't bias means
                load[load == 0] = 'nan'
                mean_load = np.nanmean(load)
                # add to dictionary
                wasielewski24_load_dict[wwtp] = mean_load

        # Mohamedali et al., 2020
        for wwtp in mohamedali2020_WWTP_ds['name'].values:
                # get mean loading of each WWTP
                load = mohamedali2020_WWTP_ds.sel(source=(mohamedali2020_WWTP_ds['name']==wwtp))['load [kg/d]'].values
                # replace all zeros with nans, so that they don't bias means
                load[load == 0] = 'nan'
                mean_load = np.nanmean(load)
                # add to dictionary
                mohamedali2020_load_dict[wwtp] = mean_load

        # sort dictionary from highest to lowest dischargers
        wasielewski24_load_dict = dict(sorted(wasielewski24_load_dict.items(), key=lambda item: item[1], reverse=True))
        mohamedali2020_load_dict = dict(sorted(mohamedali2020_load_dict.items(), key=lambda item: item[1], reverse=True))

        # plot ranked bar chart for each set of data
        fig, ax = plt.subplots(2,1,figsize = (8,8))
        # plot mohamedali et al., 2020 data
        n=10
        canadian = 'crimson'
        both_datasets = 'royalblue'
        colors = [canadian,both_datasets,both_datasets,canadian,
                  both_datasets,canadian,canadian,both_datasets,
                  both_datasets,canadian]
        ax[0].bar(range(n), list(mohamedali2020_load_dict.values())[0:n],color=colors)
        ax[0].set_xticks(range(n), list(mohamedali2020_load_dict.keys())[0:n],
                         rotation=30,ha='right')
        # plot wasielewski et al., 2024 data
        both_datasets = 'royalblue'
        was_only = 'orange'
        colors = [both_datasets,both_datasets,both_datasets,both_datasets,
                  both_datasets,was_only,both_datasets,both_datasets,
                  both_datasets,both_datasets]
        ax[1].bar(range(n), list(wasielewski24_load_dict.values())[0:n],color=colors)
        ax[1].set_xticks(range(n), list(wasielewski24_load_dict.keys())[0:n],
                         rotation=30,ha='right')
        # add color legends
        ax[0].text(0.9,0.8,'Included in both datasets',fontsize=10,fontweight='bold',
                   transform=ax[0].transAxes, horizontalalignment='right',
                   color='royalblue')
        ax[0].text(0.9,0.7,'Canadian (Mohamedali et al., 2020 only)',fontsize=10,fontweight='bold',
                   transform=ax[0].transAxes, horizontalalignment='right',
                   color='crimson')
        ax[0].text(0.9,0.6,'Wasiewlewski et al., 2024 only',fontsize=10,fontweight='bold',
                   transform=ax[0].transAxes, horizontalalignment='right',
                   color='orange')
        # format figure
        ax[0].set_ylim([0,16000])
        ax[1].set_ylim([0,16000])
        ax[0].set_ylabel('Average DIN load [kg/d]',fontsize=12)
        ax[1].set_ylabel('Average DIN load [kg/d]',fontsize=12)
        ax[0].text(0.03,0.9,'(a) Mohamedali et al., 2020',fontsize=12,fontweight='bold',
                   transform=ax[0].transAxes, horizontalalignment='left')
        ax[1].text(0.03,0.9,'(b) Wasielewski et al., 2024',fontsize=12,fontweight='bold',
                   transform=ax[1].transAxes, horizontalalignment='left')
        plt.suptitle('{} largest WWTP dischargers'.format(str(n)),fontweight='bold',fontsize=14)
        plt.tight_layout()

# plot WWTP outfall locations on an interactive map ----------------------------------
if unique_plant_locations == True:
        fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 3]},figsize = (15,5))
        # get the grid data
        ds = xr.open_dataset('../../../../../LO_data/grids/cas7/grid.nc')
        z = -ds.h.values
        mask_rho = np.transpose(ds.mask_rho.values)
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        X = lon[0,:] # grid cell X values
        Y = lat[:,0] # grid cell Y values
        plon, plat = pfun.get_plon_plat(lon,lat)
        # make a version of z with nans where masked
        zm = z.copy()
        zm[np.transpose(mask_rho) == 0] = np.nan
        zm[np.transpose(mask_rho) != 0] = -1
        # add land and water mask to both subplots
        ax[0].pcolormesh(plon, plat, zm, vmin=-20, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # open dataset of WWTP names for both datasets
        same_wwtp_names = Ldir['data'] / trapsD / 'wwtp_names.xlsx'
        same_WWTPs_df = pd.read_excel(same_wwtp_names)

        # get list of WWTPs that are in both datasets
        both_WWTPs = same_WWTPs_df['Mohamedali et al., 2020'][same_WWTPs_df['Wasielewski et al., 2024'].notna()].tolist()
        # get list of WWTPs that are only in Mohamedali et al., 2020
        moh20_WWTPs = [wwtp for wwtp in mohamedali2020_WWTP_ds['name'].values if wwtp not in both_WWTPs]
        # get list of WWTPs that are only in Wasielewski et al., 2024
        was24_WWTPs = [wwtp for wwtp in WWTP_df['FAC_NAME'] if wwtp not in same_WWTPs_df['Wasielewski et al., 2024'].values]

        # add WWTPs that are in both datasets
        lats = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(both_WWTPs))['lat'].values
        lons = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(both_WWTPs))['lon'].values
        ax[0].scatter(lons,lats, color='royalblue',s=30,alpha=0.8,
                zorder=5,label='Present in Both Datasets')
        
        # Mohamedali et al. (2020) only
        lats = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(moh20_WWTPs))['lat'].values
        lons = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(moh20_WWTPs))['lon'].values
        ax[0].scatter(lons,lats, color='crimson',s=30,alpha=0.8, marker='d',
                zorder=5,label='Mohamedali et al. (2020) only')
        
        # Wasielewski et al. (2024) only
        lats = WWTP_df.loc[WWTP_df['FAC_NAME'].isin(was24_WWTPs),'FAC_LAT'].values
        lons = WWTP_df.loc[WWTP_df['FAC_NAME'].isin(was24_WWTPs),'FAC_LON'].values
        ax[0].scatter(lons,lats, color='orange',s=30,alpha=0.8, marker='s',
                zorder=5,label='Wasielewski et al. (2024) only')
        
        # format figure
        pfun.dar(ax[0])
        pfun.add_coast(ax[0],color='paleturquoise')
        ax[0].set_title('WWTP locations',fontsize=14)
        ax[0].set_ylim([46.8,49.4])
        ax[0].set_xlim([-125.2,-121])
        ax[0].legend(loc='lower left',fontsize=10)
        
        # -----------------------------------------------------------
        # plot nutrient loads (monthly)
        t = pd.date_range(start='2005-01-01', end='2005-12-01', freq='MS')

        # Original Mohamedali et al. (2020)
        # get sum of load from all sources
        all_WWTP_loads = mohamedali2020_WWTP_ds['load [kg/d]'].sum(dim='source')
        # then, get the mean annual loading profile with monthly resolution (still, sum of all plants)
        currentLO = all_WWTP_loads.groupby('date.month').mean(dim='date')
        # plot time series of loads
        ax[1].plot(t,currentLO, color='crimson',linewidth=4,alpha=0.7,
                   label='[A] = [B+C] Current LiveOcean (all WWTPs in Moh20)')

        # Mohamedali et al. (2020) only
        # get all loading data for WWTPs in this dataset only
        WWTP_load_ds = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(moh20_WWTPs))
        # get sum of load from all sources
        all_WWTP_loads = WWTP_load_ds['load [kg/d]'].sum(dim='source')
        # then, get the mean annual loading profile with monthly resolution (still, sum of all plants)
        Moh20_only = all_WWTP_loads.groupby('date.month').mean(dim='date')
        # plot time series of loads
        ax[1].plot(t,Moh20_only, color='crimson',linewidth=3,alpha=0.7,
                   linestyle=':', label='[B] WWTPs only in Moh20')
        
        # Both datasets, but Mohamedali et al. (2020) version
        # get all loading data for WWTPs in this dataset only
        WWTP_load_ds = mohamedali2020_WWTP_ds.sel(source=mohamedali2020_WWTP_ds['name'].isin(both_WWTPs))
        # get sum of load from all sources
        all_WWTP_loads = WWTP_load_ds['load [kg/d]'].sum(dim='source')
        # then, get the mean annual loading profile with monthly resolution (still, sum of all plants)
        both_Moh20ver = all_WWTP_loads.groupby('date.month').mean(dim='date')
        # plot time series of loads
        ax[1].plot(t,both_Moh20ver, color='royalblue',linewidth=3,alpha=0.7,
                   linestyle='--', label='[C] WWTPs in both datasets (Moh20 loads)')
        
        # Both datasets, but Wasielewski et al. (2024) version
        # get all loading data for WWTPs in this dataset only
        both_wwtps_was24names = [wwtp for wwtp in WWTP_df['FAC_NAME'] if wwtp in same_WWTPs_df['Wasielewski et al., 2024'].values]
        IDlist = WWTP_df.loc[WWTP_df['FAC_NAME'].isin(both_wwtps_was24names),'FAC_ID'].values
        WWTP_load_df = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'].isin(IDlist)]
        # first, sum all of the individual plants together to get one 2005 - 2020 timeseries (which is the sum of all plants)
        all_WWWTP_loads = WWTP_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        # then, get the mean annual loading profile (still, sum of all plants)
        annualmean_all_WWTP_loads = all_WWWTP_loads.groupby(['MONTH']).mean()
        # plot time series of loads
        month_days = [1/31, 1/28.25, 1/31, 1/30,
                      1/31, 1/30,    1/31, 1/31,
                      1/30, 1/31,    1/30, 1/31]
        # plot time series of loads
        both_Was24ver = annualmean_all_WWTP_loads*month_days
        ax[1].plot(t,both_Was24ver, color='royalblue',linewidth=3,alpha=0.7,
                   linestyle='-', label='[D] WWTPs in both datasets (Was24 loads)')
        
        # WWTPs in Wasielewski et al. (2024) only
        # get all loading data for WWTPs in this dataset only
        IDlist = WWTP_df.loc[WWTP_df['FAC_NAME'].isin(was24_WWTPs),'FAC_ID'].values
        WWTP_load_df = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'].isin(IDlist)]
        # first, sum all of the individual plants together to get one 2005 - 2020 timeseries (which is the sum of all plants)
        all_WWWTP_loads = WWTP_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        # then, get the mean annual loading profile (still, sum of all plants)
        annualmean_all_WWTP_loads = all_WWWTP_loads.groupby(['MONTH']).mean()
        # plot time series of loads
        month_days = [1/31, 1/28.25, 1/31, 1/30,
                      1/31, 1/30,    1/31, 1/31,
                      1/30, 1/31,    1/30, 1/31]
        # plot time series of loads
        Was24_only = annualmean_all_WWTP_loads*month_days
        ax[1].plot(t,Was24_only, color='orange',linewidth=3,alpha=0.7,
                   linestyle=':', label='[E] WWTPs only in Was24')
        
        # All WWTPs in Wasielewski et al. (2024)
        # get all loading data for WWTPs in this dataset only
        IDlist = WWTP_df['FAC_ID'].values
        WWTP_load_df = wasielewski24_WWTP_df.loc[wasielewski24_WWTP_df['FAC_ID'].isin(IDlist)]
        # first, sum all of the individual plants together to get one 2005 - 2020 timeseries (which is the sum of all plants)
        all_WWWTP_loads = WWTP_load_df.groupby(['YEAR', 'MONTH'])['TN_LOAD_KG_MO'].sum()
        # then, get the mean annual loading profile (still, sum of all plants)
        annualmean_all_WWTP_loads = all_WWWTP_loads.groupby(['MONTH']).mean()
        # plot time series of loads
        month_days = [1/31, 1/28.25, 1/31, 1/30,
                      1/31, 1/30,    1/31, 1/31,
                      1/30, 1/31,    1/30, 1/31]
        # plot time series of loads
        all_Was24 = annualmean_all_WWTP_loads*month_days
        ax[1].plot(t,all_Was24, color='orange',linewidth=3,alpha=0.7,
                   label='[F] = [D+E] all WWTPs in Was24')
        
        # Proposed dataset
        ax[1].plot(t,both_Was24ver+Moh20_only, color='black',linewidth=4,alpha=0.7,
                   label='[G] = [B+D] Proposed')

        # format figure
        ax[1].set_title('Mean Annual Nutrient Loads',fontsize=14)
        ax[1].set_ylabel('Total Nitrogen Load [kg/day]',fontsize=12)
        ax[1].set_ylim([0,90000])
        ax[1].set_xlim([t[0],t[-1]])
        ax[1].xaxis.set_major_locator(mdates.MonthLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax[1].legend(loc='best',ncol=2)

        plt.tight_layout()
        plt.show()

# TODO: verify the total N load variable is the same as if I manually calculate it (it's not, and I need to correct the above figures)
# TODO: save figures in this section
# TODO: incorporate new_wwtp_loading.py plotting
# TODO: check loads in Hood Canal from WWTPs and fish hatcheries