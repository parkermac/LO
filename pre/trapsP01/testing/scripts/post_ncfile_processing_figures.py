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
#                           Create output directory                             #
#################################################################################

plt.close('all')

out_dir = Ldir['LOo'] / 'point_source_integration' / 'final_WWTP_loads'
Lfun.make_dir(out_dir)

print('\n')

#################################################################################
#                   Wasielewski et al., 2024 loading plots                      #
#################################################################################

plot_individual_WWTP_loads = False
plot_facility_type_load_comparison = False
plot_hoodcanal_load_comparison = True

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent.parent.parent
with open(this_dir /'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# location of point source data to process
data_dir = Ldir['data'] / trapsD / 'processed_data'
was24_data_fn = data_dir / 'wwtp_data_wasielewski_etal_2024.nc'

# get point source data
was24_data = xr.open_dataset(was24_data_fn)

#################################################################################
#                    Mohamedali et al., 2020 loading plots                      #
#################################################################################

# location of point source data to process
moh20_data_fn = data_dir / 'wwtp_data_mohamedali_etal_2020.nc'
# get point source data
moh20_data = xr.open_dataset(moh20_data_fn)

#################################################################################
#                             Plot proposed loads                               #
#################################################################################

        
# loop through all WWTPs in Mohamedali et al. (2020):
for i,wwtp in enumerate(moh20_data['name'].values):
        
        # initialize figure
        fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 3]},figsize = (14,5))
        
        # plot wwtp location on map -------------------------------------------
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

        # plot wwtp location
        lat = moh20_data.sel(source=moh20_data['name']==wwtp)['lat']
        lon = moh20_data.sel(source=moh20_data['name']==wwtp)['lon']
        ax[0].scatter(lon,lat,color='navy',s=100,marker='*',zorder=3)

        # format figure
        pfun.dar(ax[0])
        pfun.add_coast(ax[0],color='paleturquoise')
        ax[0].set_title('{} location'.format(wwtp),fontsize=14)
        ax[0].set_ylim([46.8,49.4])
        ax[0].set_xlim([-125,-121.5])

        # add loading profiles --------------------------------------------------------

        # Subsample to one value per month
        moh20_data['load [kg/d]'] = moh20_data['flow'] * (moh20_data['NO3'] + moh20_data['NH4']) * 14.01 * 60 * 60 * 24 / 1000 / 1000
        moh20_raw = moh20_data.sel(source=(moh20_data['name']==
                                wwtp))['load [kg/d]'].resample(date='1MS').first()[0]
        
        
        # get climatology if using old data (i.e. WWTP not available in was24):
        # get flow [m3/s]
        # fn_oldWWTP_flow = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_flow.p'
        fn_oldWWTP_flow = Ldir['LOo'] / 'pre' / 'trapsP01' / 'moh20_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_flow.p'
        df_oldWWTP_flow = pd.read_pickle(fn_oldWWTP_flow)
        # Get WWTP flow
        old_flow = df_oldWWTP_flow[wwtp].values # [m3/s]

        # get NO3 and NH4 [mmol/m3]
        # fn_oldWWTP_NO3 = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_NO3.p'
        fn_oldWWTP_NO3 = Ldir['LOo'] / 'pre' / 'trapsP01' / 'moh20_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_NO3.p'
        df_oldWWTP_NO3 = pd.read_pickle(fn_oldWWTP_NO3)
        # fn_oldWWTP_NH4 = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_NH4.p'
        fn_oldWWTP_NH4 = Ldir['LOo'] / 'pre' / 'trapsP01' / 'moh20_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_NH4.p'
        df_oldWWTP_NH4 = pd.read_pickle(fn_oldWWTP_NH4)
        # Get WWTP NO3 and NH4
        old_NO3 = df_oldWWTP_NO3[wwtp].values # [mmol/m3]
        old_NH4 = df_oldWWTP_NH4[wwtp].values # [mmol/m3]
        # sum nutrients
        old_nutrients = old_NO3 + old_NH4 # [mmol/m3]

        # get total nutrient loading [kg/day]
        old_nutrient_load_mmol_s = old_flow * old_nutrients # [m3/s * mmol/m3 = mmol/s]
        # convert to kg/day: [1g = 14.01/1000^2mmol] [1day = 60*60*24sec]
        old_nutrient_load_kg_day = old_nutrient_load_mmol_s * 14.01 * 60 * 60 * 24 / 1000 / 1000

        # get one value every month
        old_nutrient_load_subsampled = old_nutrient_load_kg_day[15::30]

        # repeat for 15 years
        moh20_climatology = np.tile(old_nutrient_load_subsampled,27)

        # get time array
        t_old = pd.date_range(start='1999-01-01', end='2017-7-31', freq='MS').to_list()
        t_all = pd.date_range(start='1999-01-01',end='2025-12-31', freq='MS').to_list()


        # plot old WWTP climatologies
        ax[1].plot(t_all,moh20_climatology, label='Climatology (based on Mohamedali et al., 2020)',
                linewidth=3,color='k', alpha=0.4)
                
        # plot Mohamedali et al. (2020) raw data
        ax[1].plot(t_old,moh20_raw, label='Mohamedali et al. (2020) raw',
                linewidth=1.5, linestyle='--', color='purple', alpha=0.8)

        # format figure
        ax[1].set_title('{} nutrient load comparison'.format(wwtp),
                fontsize=14,fontweight='bold')
        ax[1].set_ylabel('DIN nutrient load [kg/day]',fontsize=12)
        ax[1].legend(loc='best', fontsize=12)

        ax[1].set_ylim([0,np.nanmax(moh20_raw)*1.1])
        ax[1].set_xlim([pd.to_datetime('1999-01-01'),pd.to_datetime('2025-12-31')])
        ax[1].xaxis.set_major_locator(mdates.YearLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax[1].tick_params(axis='x', labelrotation=30)
        ax[1].grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')

        # save figure ---------------------------------------------------------------
        plt.tight_layout()
        plt.savefig(out_dir / (wwtp + '.png'))

        plt.close()


# loop through all WWTPs in Wasielewski et al. (2024):
for i,wwtp_ID in enumerate(was24_data['ID'].values):
        
        # initialize figure
        fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 3]},figsize = (14,5))
        
        # plot wwtp location on map -------------------------------------------
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

        # plot wwtp location
        lat = was24_data.sel(source=was24_data['ID']==wwtp_ID)['lat']
        lon = was24_data.sel(source=was24_data['ID']==wwtp_ID)['lon']
        ax[0].scatter(lon,lat,color='navy',s=100,marker='*',zorder=3)

        # format figure
        pfun.dar(ax[0])
        pfun.add_coast(ax[0],color='paleturquoise')
        wwtp_name = was24_data.sel(source=was24_data['ID']==wwtp_ID)['name'].values[0]
        ax[0].set_title('{} location'.format(wwtp_name),fontsize=14)
        ax[0].set_ylim([46.8,49.4])
        ax[0].set_xlim([-125,-121.5])

        # get climatologies -----------------------------------------------------------

        # get flow [m3/s]
        fn_oldWWTP_flow = Ldir['LOo'] / 'pre' / 'trapsP01' / 'was24_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_flow.p'
        df_oldWWTP_flow = pd.read_pickle(fn_oldWWTP_flow)
        # Get WWTP flow
        old_flow = df_oldWWTP_flow[wwtp_name].values # [m3/s]

        # get NO3 and NH4 [mmol/m3]
        fn_oldWWTP_NO3 = Ldir['LOo'] / 'pre' / 'trapsP01' / 'was24_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_NO3.p'
        df_oldWWTP_NO3 = pd.read_pickle(fn_oldWWTP_NO3)
        fn_oldWWTP_NH4 = Ldir['LOo'] / 'pre' / 'trapsP01' / 'was24_wwtps' / 'lo_base' / 'Data_historical' / 'CLIM_NH4.p'
        df_oldWWTP_NH4 = pd.read_pickle(fn_oldWWTP_NH4)
        # Get WWTP NO3 and NH4
        old_NO3 = df_oldWWTP_NO3[wwtp_name].values # [mmol/m3]
        old_NH4 = df_oldWWTP_NH4[wwtp_name].values # [mmol/m3]
        # sum nutrients
        old_nutrients = old_NO3 + old_NH4 # [mmol/m3]

        # get total nutrient loading [kg/day]
        old_nutrient_load_mmol_s = old_flow * old_nutrients # [m3/s * mmol/m3 = mmol/s]
        # convert to kg/day: [1g = 14.01/1000^2mmol] [1day = 60*60*24sec]
        old_nutrient_load_kg_day = old_nutrient_load_mmol_s * 14.01 * 60 * 60 * 24 / 1000 / 1000

        # get one value every month
        old_nutrient_load_subsampled = old_nutrient_load_kg_day[15::30]

        # repeat for 15 years
        was24_climatology = np.tile(old_nutrient_load_subsampled,5)

        # add loading profiles --------------------------------------------------------

        # Subsample to one value per month
        was24_data['load [kg/d]'] = was24_data['flow'] * (was24_data['NO3'] + was24_data['NH4']) * 14.01 * 60 * 60 * 24 / 1000 / 1000
        was24_raw = was24_data.sel(source=(was24_data['ID']==
                                wwtp_ID))['load [kg/d]'].resample(date='1MS').first()[0]

        # get time array
        t_new = pd.date_range(start='2005-01-01', end='2020-12-31', freq='MS').to_list()
        t_wasclim = pd.date_range(start='2021-01-01',end='2025-12-31', freq='MS').to_list()

                
        # get Wasielewski et al. (2024) data
        ax[1].plot(t_new,was24_raw, label='Wasielewski et al. (2024) raw',
                                linewidth=2.5,color='hotpink', alpha=0.8) 

        # plot WWTP climatologies
        ax[1].plot(t_wasclim,was24_climatology, label='Climatology (based on Wasielewski et al., 2024)',
                linewidth=3,color='k', alpha=0.4)   

        # format figure
        ax[1].set_title('{} nutrient load comparison'.format(wwtp_name),
                fontsize=14,fontweight='bold')
        ax[1].set_ylabel('DIN nutrient load [kg/day]',fontsize=12)
        ax[1].legend(loc='best', fontsize=12)

        ax[1].set_ylim([0,np.nanmax(was24_raw)*1.1])
        ax[1].set_xlim([pd.to_datetime('1999-01-01'),pd.to_datetime('2025-12-31')])
        ax[1].xaxis.set_major_locator(mdates.YearLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax[1].tick_params(axis='x', labelrotation=30)
        ax[1].grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')

        # save figure ---------------------------------------------------------------
        plt.tight_layout()
        plt.savefig(out_dir / (wwtp_name + '.png'))
        plt.close()