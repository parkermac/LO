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

plot_individual_WWTP_loads = True
plot_facility_type_load_comparison = True

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent.parent.parent
with open(this_dir /'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# location of point source data to process
pointsource_dir = Ldir['data'] / trapsD / 'raw_data' / 'point_sources'
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
        ax[0].set_ylim([46.8,49])
        ax[0].set_xlim([-124.6,-121])
        ax[0].legend(loc='lower left',fontsize=14)

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
        ax[1].plot(t,annualmean_all_WWTP_loads, color='deeppink',linewidth=3,alpha=0.7)
        ax[1].plot(t,annualmean_all_hatchery_loads, color='black',linewidth=3,alpha=0.7)
        ax[1].plot(t,annualmean_all_industrial_loads, color='deepskyblue',linewidth=3,alpha=0.7)
        # format figure
        ax[1].set_title('Mean Annual Nutrient Loads',fontsize=14)
        ax[1].set_ylabel('Total Nitrogen Load\n[kg/month]',fontsize=12)
        ax[1].set_yscale('log')
        ax[1].set_ylim([1,1e7])
        ax[1].set_xlim([t[0],t[-1]])
        ax[1].xaxis.set_major_locator(mdates.MonthLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

        # format and save figure
        plt.suptitle('Wasielewski et al. (2024)',fontsize=14,fontweight='bold')
        plt.tight_layout
        plt.show()
        plt.savefig(out_dir / 'Wasielewski_etal2024' / 'all_loads_comparison.png')


#################################################################################
#                   Wasielewski et al., 2024 loading plots                      #
#################################################################################