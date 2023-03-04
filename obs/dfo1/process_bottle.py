"""
Code to process DFO bottle data from the NetCDF version.

Performance:
"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time
import xarray as xr

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()
in_dir = Ldir['data'] / 'obs' / 'dfo1'

# get the dataset
fn = in_dir / 'IOS_BOT_Profiles.nc'
print(fn.name)
ds = xr.open_dataset(fn,decode_times=False)

testing = False

# initial exploration
if False:
    for vn in ds.data_vars:
        try:
            units = ds[vn].units
        except:
            units = ''
        try:
            long_name = ds[vn].long_name
        except:
            long_name = ''
        print('%s: %s [%s]' % (vn,long_name,units))
"""
RESULT (* = variables we will use):
CNDCST01: Sea Water Electrical Conductivity [S/m]
*CPHLFLP1: Concentration of chlorophyll-a per unit volume of the water body [mg/m^3]
*DOXMZZ01: Oxygen concentration [umol/kg]
DOXYZZ01: Oxygen concentration [mL/L]
*NTRZAAZ1: Mole Concentration of Nitrate and Nitrite in Sea Water [umol/L]
*PHOSAAZ1: Mole Concentration of Phosphate in Sea Water [umol/L]
PRESPR01: Pressure [decibar]
PSALBST01: Sea Water Practical Salinity [PSS-78]
PSALST01: Sea Water Practical Salinity [PSS-78]
PSALST02: Sea Water Practical Salinity [PSS-78]
*SLCAAAZ1: Mole Concentration of Silicate in Sea Water [umol/L]
SSALBST01: Sea Water Salinity [PPT]
SSALST01: Sea Water Salinity [PPT]
TEMPRTN1: Sea Water Temperature [degC]
TEMPS601: Sea Water Temperature [degC]
TEMPS602: Sea Water Temperature [degC]
TEMPS901: Sea Water Temperature [degC]
TEMPS902: Sea Water Temperature [degC]
TEMPST01: Sea Water Temperature [degC]
agency:  []
country:  []
*depth: Depth [m]
*event_number:  []
filename:  []
geographic_area:  []
instrument_model:  []
instrument_serial_number:  []
instrument_type:  []
*latitude: Latitude [degrees_north]
*longitude: Longitude [degrees_east]
mission_id:  []
platform:  []
profile: Profile ID []
project:  []
scientist:  []
*sea_water_practical_salinity: Sea Water Practical Salinity [PSS-78]
sea_water_pressure: Pressure [dbar]
*sea_water_temperature: Sea Water Temperature [degC]
*time: Time []
"""
col_dict = {'CPHLFLP1':'Chl (mg m-3)','DOXMZZ01': 'DO (umol kg-1)','NTRZAAZ1':'NO3 (uM)',
    'PHOSAAZ1':'PO4 (uM)','SLCAAAZ1':'SiO4 (uM)',
    'depth':'depth','event_number':'cid','latitude':'lat','longitude':'lon',
    'sea_water_practical_salinity':'SP','sea_water_temperature':'TI','time':'tsec',
    # 'sea_water_pressure':'pres',
    }

if testing:
    year_list = [2017]
else:
    year_list = list(range(1930, datetime.now().year + 1))
out_dir = Ldir['LOo'] / 'obs' / 'dfo1' / 'bottle'
Lfun.make_dir(out_dir)

tt0 = time()
dt_ref = datetime(1970,1,1)
for year in year_list:
    out_fn = out_dir / (str(year) + '.p')
    info_out_fn = out_dir / ('info_' + str(year) + '.p')

    # limit time range
    dt0 = datetime(year,1,1)
    dt1 = datetime(year+1,1,1)
    tsec0 = (dt0 - dt_ref).total_seconds()
    tsec1 = (dt1 - dt_ref).total_seconds()
    mask = (ds.time.values>=tsec0) & (ds.time.values<tsec1)
    
    if mask.sum() > 0:
    
        # read selected variables into a DataFrame
        df = pd.DataFrame()
        for vn in col_dict.keys():
            df[vn] = ds[vn][mask]
        df = df.rename(col_dict, axis=1)
    
        # there are a few fields that had values like 1e37
        for vn in ['depth','TI','SP']:
            df.loc[df[vn]>1e5,vn] = np.nan
        # and a few negative depths and salinities
        df[df['depth']<=0] = np.nan
        df[df['SP']<=0] = np.nan
        
        # Exploration: try to get z two ways
        # df.loc[df['press']>1e5,'press'] = np.nan
        # df[df['pres']<=0] = np.nan
        # z1 = -df.depth.to_numpy()
        # z2 = gsw.z_from_p(df.pres.to_numpy(), df.lat)
        # nz1 = (~np.isnan(z1)).sum()
        # nz2=(~np.isnan(z2)).sum()
        # print('%d: z from depth = %d, z from pres = %d (diff=%d)' % (year,nz1,nz2,nz1-nz2))
        # RESULT: mostly the same number, occasionally the one from depth has a few more.
        # We will choose to just use z from depth
        
        # add z and drop rows missing any fundamental variables
        df['z'] = -df.depth
        df = df.dropna(subset=['z','lon','lat'])
            
        # convert to standard units
        p = gsw.p_from_z(df.z.to_numpy(), df.lat) # no clue why .to_numpy() is needed
        df['SA'] = gsw.SA_from_SP(df.SP, p, df.lon, df.lat)
        df['CT'] = gsw.CT_from_t(df.SA, df.TI, p)
        rho = gsw.rho(df.SA, df.CT, p)
        df['DO (uM)'] = df['DO (umol kg-1)'].to_numpy() * rho / 1000
        # NOTE: convenient ways to look for nans or out of range values are
        # df.isna().sum() and df.max(), df.min()
        tsec = df.tsec.to_numpy()
        time_list = []
        for item in tsec:
            time_list.append(dt_ref + timedelta(seconds=item))
        df['time'] = time_list

        # clean up columns
        cols = ['cid', 'lon', 'lat', 'time', 'z','SA', 'CT',
            'DO (uM)', 'NO3 (uM)', 'NH4 (uM)', 'PO4 (uM)' 'SiO4 (uM)', 'Chl (mg m-3)']
        this_cols = [item for item in cols if item in df.columns]
        df = df[this_cols]

        df['name'] = None
        df['cruise'] = None

        # Renumber cid to be increasing from zero in steps of one.
        df = obs_functions.renumber_cid(df)

        if len(df) > 0:
            # Save the data
            df.to_pickle(out_fn)
            info_df = obs_functions.make_info_df(df)
            info_df.to_pickle(info_out_fn)
            print('%d: %d casts' % (year,int(info_df.index.max())))
        else:
            print('%d: 0 casts' % (year))
    else:
        print('%d: 0 casts' % (year))


print('Elapsed time for bottle = %0.1f sec' % (time()-tt0))

