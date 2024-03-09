"""
- Code to process the Line P complied bottle data to pickle files for year 1993-1998
The data range from year 1990-2019.
However,
the complied dataset doesn't have consistent columns for temp, salinity, even in the same year (e.g.,1994)

During this period "TMP_REVERSING_DEG_C", "CTDTMP_ITS90_DEG_C", "SALINITY_PSS78", "CTDSAL_PSS78" are cross-used

Some notes:

- year 1993 doesn't have P4, but has P12
- year 1995 doesn't have P4, nor P12
- year 1998 doesn't have P4, nor P12

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'LineP'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source / otype

testing = False

if testing:
    year_list = [1998]
else:
    year_list = range(1993, 1998+1)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save. I looked
# at units_dict (created below) to be sure about the units.
v_dict = {
    'time':'time',
    'EXPOCODE':'',
    'CRUISE_ID':'cruise',
    'STATION_ID':'name',
    'EVENT_NO':'cid',
    'NISKIN_NO':'',
    'YEAR_UTC':'',
    'MONTH_UTC':'',
    'DAY_UTC':'',
    'TIME_UTC':'',
    'YEARDAY_UTC':'',
    'LONGITUDE_DEC':'lon',
    'LATITUDE_DEC':'lat',
    'CTDPRS_DBAR':'P (dbar)',
    'CTDTMP_ITS90_DEG_C':'IT1',
    'TMP_REVERSING_DEG_C':'IT2',
    'CTDSAL_PSS78':'SP1',
    'CTDSAL_FLAG_W':'',
    'SALINITY_PSS78':'SP2',
    'SALINITY_FLAG_W':'',
    'OXYGEN_UMOL_KG':'DO (umol/kg)',
    'OXYGEN_FLAG_W':'',
    'DIC_UMOL_KG':'DIC (umol/kg)',
    'DIC_FLAG_W':'',
    'TA_UMOL_KG':'TA (umol/kg)',
    'TA_FLAG_W':'',
    'NITRATE_NITRITE_UMOL_KG':'NO3 (umol/kg)',
    'NITRATE_NITRITE_FLAG_W':'',
    'SILICATE_UMOL_KG':'SiO4 (umol/kg)',
    'SILICATE_FLAG_W':'',
    'PHOSPHATE_UMOL_KG':'PO4 (umol/kg)',
    'PHOSPHATE_FLAG_W':'',
    'original_filename':'',
}

load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)

    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')

    if year in year_list:
        in_fn = in_dir0 / 'LineP_for_Data_Synthesis_1990-2019_v1.csv'
        """
        there are two rows have bad data, which is row 700, and 724, omit these 
        """
        df0 = pd.read_csv(in_fn, low_memory=False,
                          parse_dates={'time':['YEAR_UTC','MONTH_UTC','DAY_UTC','TIME_UTC']},
                          usecols= list(range(33)))
        df0 = df0.drop([699, 723])
        df0 = df0.reset_index(drop=True)
        df0['time'] = pd.to_datetime(df0['time'], format='ISO8601')
        load_data = False # only load the first time
        # select one year
    t = pd.DatetimeIndex(df0.time)
    df1 = df0.loc[t.year==year,:].copy()

    # select and rename variables
    df = pd.DataFrame()
    for v in df1.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df[v_dict[v]] = df1[v]

    # missing data is -999
    df[df==-999] = np.nan

    # a little more cleaning up
    df = df.dropna(axis=0, how='all') # drop rows with no good data
    df = df[df.time.notna()] # drop rows with bad time
    df = df.reset_index(drop=True)

    """
    Line P goes way beyond the model domain, now limit the geographical region to the model domain, 
    most of the time only P4 is with the domain,
    but P12 is at the edge of the model domain, -130.67 
    So, include this station because we can use it to test the boundary condition
    """
    df = df[(df.lon>-131) & (df.lon<-122) & (df.lat>42) & (df.lat<52)]

    # This dataset doesn't have unique numbers for each cast, if same stations,
    # force lat and lon to be consistent throughout the given year
    for name in df.name.unique():
        df.loc[df.name==name,'lon'] = df[df.name==name].lon.values[0]
        df.loc[df.name==name,'lat'] = df[df.name==name].lat.values[0]

    # Each cast is associated with a different time (hh:mm), now assign arbitrary cid based on this
    unique = df.time.unique()
    ind = {time: cid for time, cid in zip(unique, range(len(unique)))}
    df['cid'] = df['time'].map(ind)

    # Next make derived quantities and do unit conversions
    # (1) Create CT, SA, and z
    # - pull out variables
    SP2 = df.SP2.to_numpy()
    IT2 = df.IT2.to_numpy()
    SP1 = df.SP1.to_numpy()
    IT1 = df.IT1.to_numpy()
    # combine
    SP = np.where(np.isnan(SP1), SP2, SP1)
    IT = np.where(np.isnan(IT1), IT2, IT1)
    p = df['P (dbar)'].to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()

    # - do the conversions
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    z = gsw.z_from_p(p, lat)
    # - add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    df['z'] = z
    rho = gsw.rho(SA, CT, p)

    # (2) units
    for vn in ['DO', 'NO3', 'PO4', 'SIO4', 'TA', 'DIC']:
        if (vn+' (umol/kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol/kg)']

    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
            'CT', 'SA', 'DO (uM)',
            'NO3 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
            'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    print(' - processed %d casts' % (len(df.cid.unique())))

    # Renumber cid to be increasing from zero in steps of one.
    df = obs_functions.renumber_cid(df)

    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)
