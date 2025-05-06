"""
Script to plot and compare new WWTP data and my older climatologies.
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
import pandas as pd
import argparse
import datetime as dt
import matplotlib.dates as mdates
import os
from pathlib import Path

#################################################################################
#                             Facility information                              #
#################################################################################

plt.close('all')
# where to put output figures
out_dir = 'testing/loading_figures'
Lfun.make_dir(out_dir)

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent
with open(this_dir / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# read csv files
pointsource_dir = Ldir['data'] / trapsD / 'raw_data' / 'point_sources'
pointsource_meta = pointsource_dir / 'fac_attributes.csv'
pointsource_loads = pointsource_dir / 'nutrient_loads.csv'
df_fac_attr = pd.read_csv(pointsource_meta)
df_nut_loads_new = pd.read_csv(pointsource_loads)

facilities_all = df_fac_attr['FAC_ID'].values.tolist()
facilities = []


# remove fish hatcheries and industrial facilities
for fac_ID in facilities_all:
        # get facility type
        fac_type = df_fac_attr.loc[df_fac_attr['FAC_ID'] == fac_ID, 'FAC_TYPE'].values[0]
        # remove anything that is not a WWTP
        if 'sic_4952' in fac_type:
                facilities = facilities + [fac_ID]

# remove WWTPs that have incomplete time series,
# or are irrelevant to LiveOcean (initialized at end of 2012)
facilities.remove('WA0020567-001') # OAK HARBOR STP which ends mid-way through 2010


for fac_ID in facilities:

        #################################################################################
        #                               Get new WWTP data                               #
        #################################################################################

        # get WWTP ID
        fac_ID = fac_ID

        # get name
        fac_name = df_fac_attr.loc[df_fac_attr['FAC_ID'] == fac_ID, 'FAC_NAME'].values[0]

        # get flow corresponding to desired WWTP
        # replace any '.' with NaN
        no_empties_flow = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'FLOW_MGD'].replace('.', np.nan).astype(float)
        new_flow_millGall_day = no_empties_flow.values
        # convert to m3/s
        new_flow_m3_s = new_flow_millGall_day * 0.0438126364

        # get nutrient data [mg/L]
        no_empties_NO3 = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'NO2NO3N_MG_L'].replace('.', np.nan).astype(float)
        no_empties_NH4 = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'NH4N_MG_L'].replace('.', np.nan).astype(float)
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


        #################################################################################
        #                                     Plot                                      #
        #################################################################################

        # get time array
        # t_new = pd.date_range(start='2005-01-01', end='2020-12-31', freq='MS').to_list()
        month = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'MONTH'].replace('.', np.nan).astype(int)
        year = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'YEAR'].replace('.', np.nan).astype(int)
        # create time array
        t_new = np.array([dt.datetime(year, month, 1) for year, month in zip(year, month)], dtype='datetime64[M]') 

        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (13,7))

        # plot new data
        ax.plot(t_new,new_nutrient_load_15years,marker='o',
                linewidth=1,color='hotpink', alpha=0.8)

        # format figure
        ax.set_title('{} nutrient load comparison'.format(fac_name),
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

        plt.savefig(out_dir + '/' + fac_ID +'.png')