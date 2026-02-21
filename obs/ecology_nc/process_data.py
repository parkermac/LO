"""
Code to process the ecology NetCDF ctd and bottle data.

Takes 16 seconds to run for 25 years (!).

"""

import pandas as pd
import numpy as np
import xarray as xr
import gsw
import sys
from time import time as Time

from lo_tools import Lfun, zfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology_nc'
in_dir = Ldir['data'] / 'obs' / source

testing = False

if testing:
    year_list = [2023]
else:
    year_list = range(1999,2026)
    #year_list = [2005,2006,2007,2024]

# output locations are created at the end because we do both ctd and
# bottle cases from the same input file

tt0 = Time()
for year in year_list:
    ys = str(year)
    print('\n'+ys)

    # input file
    in_fn = in_dir / ('MarineWaterProfilesAndNutrientsYear' + ys + '.nc')
    ds = xr.open_dataset(in_fn)

    # These have dimension (stations)
    sta_num = ds.station_number.values    # array with values 1 to ## {*}
    sta_lon = ds.Longitude.values
    sta_lat = ds.Latitude.values
    sta_name = [item.decode() for item in ds.Station.values]

    # These have dimension (profiles)
    pro_sta_num = ds.station_index.values # array with values 1 to ## {*}
    pro_num = ds.profile_index.values     # array with values related to profiles {**}
    pro_date_dti = pd.to_datetime(ds.FieldDate.values)

    # These have dimension (obs)
    obs_pro_num = ds.obs_index.values     # array with values related to profiles {**}
    obs_dti = pd.to_datetime(ds.UTCDatetime.values)
    z = -ds.Depth.values

    # populate a DataFrame "df" with the data from the NetCDF file
    df = pd.DataFrame()
    df['cid'] = obs_pro_num
    df['z'] = z

    for N in range(0,len(pro_num)):
        # Fill in fields that we want to have one value for all
        # rows associated with a profile.
        this_pro_num = pro_num[N]
        this_sta_num = pro_sta_num[N]
        # get other info about this cast
        sta_ind = np.argwhere(sta_num==this_sta_num)
        sta_ind = sta_ind[0][0]
        this_sta_name = sta_name[sta_ind]
        this_lon = float(sta_lon[sta_ind])
        this_lat = float(sta_lat[sta_ind])
        # mask isolates a single cast
        mask = obs_pro_num == this_pro_num
        ti = obs_dti[mask]
        this_ti0 = ti[0] # assign the time from the start of the cast
        # add things to df
        df.loc[mask,'lon'] = this_lon
        df.loc[mask,'lat'] = this_lat
        df.loc[mask,'time'] = this_ti0
        df.loc[mask,'name'] = this_sta_name

    # Make a dictionary with names of data fields to process.
    # The keys are fields in the original NetCDF file, and the
    # values are names we use in the draft DataFrame.
    vn_dict = {'Temp':'Temp insitu',
        'Salinity':'Salinity psu',
        'DOAdjusted':'DO mg/L',
        'NO3':'NO3 (uM)',
        'NH4':'NH4 (uM)',
        'PO4':'PO4 (uM)',
        'SiOH4':'SiO4 (uM)',
        'FluorAdjusted':'Chl (mg m-3)', # Is this correct?
        }
    # add data to DataFrame
    for vn in vn_dict.keys():
        df[vn_dict[vn]] = ds[vn].values

    # derived variables
    # SA and CT
    SP = df['Salinity psu'].to_numpy()
    IT = df['Temp insitu'].to_numpy()
    z = df.z.to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    p = gsw.p_from_z(z, lat)
    # - do the conversions
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    # - add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    # Unit conversion
    df['DO (uM)'] = df['DO mg/L'] * 1000 / 32
    df['cruise'] = None
            
    # retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)','Chl (mg m-3)',
        'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
        'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    # At this point we have ALL the data, but to be consistent with
    # other items in the LO obs collection we want to split it into
    # ctd and bottle parts.
    #
    # pull out the bottles, assuming that any row with a bottle has NO3.
    bmask = df['NO3 (uM)'].notnull()
    bot_df = df.loc[bmask,:]
    # Renumber cid to be increasing from zero in steps of one.
    bot_df = obs_functions.renumber_cid(bot_df)
    print(' - processed %d bottle casts' % ( len(bot_df.cid.unique()) ))
    #
    # pull out ctd by dropping all the extra columns
    ctd_cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)', 'Chl (mg m-3)']
    this_cols = [item for item in ctd_cols if item in df.columns]
    ctd_df = df[this_cols].copy()
    # Renumber cid to be increasing from zero in steps of one.
    ctd_df = obs_functions.renumber_cid(ctd_df)
    print(' - processed %d ctd casts' % ( len(ctd_df.cid.unique()) ))
    #
    # Note that by doing separate renumber_cid() operations we may lose the connection
    # between cid numbers across the bottle and ctd version. I don't think this will
    # matter in practice.

    for otype in ['ctd','bottle']:
        if otype == 'ctd':
            dff = ctd_df
        elif otype == 'bottle':
            dff = bot_df
        if (len(dff) > 0):
            info_df = obs_functions.make_info_df(dff)
            if testing == False:
                # output locations
                out_dir = Ldir['LOo'] / 'obs' / source / otype
                Lfun.make_dir(out_dir)
                # name output files
                out_fn = out_dir / (ys + '.p')
                info_out_fn = out_dir / ('info_' + ys + '.p')
                # Save the data
                dff.to_pickle(out_fn)
                info_df.to_pickle(info_out_fn)
        else:
            print('No data in ' + ys + ' for ' + otype)
        
print('Total time = %d sec' % (int(Time()-tt0)))

