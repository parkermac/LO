"""
This is the main program for making the RIV forcing file.  It has been re-coded to make use
of the new pre/river output, and the alpha/river_functions.py.  It does not use the
old river class framework.

2021.11.19 This is identical to riv0 except in the calculation of Vshape.  This uses even
weights per bin which is a better approximation when using Vtranform = 2.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -test True -f riv0

run make_forcing_main.py -g cas6 -t v3 -r forecast -s continuation -d [Today's date] -test False -f riv0

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************
import os
import pandas as pd

from lo_tools import zrfun
import rivfun
from importlib import reload
reload(rivfun)

# set up the time index for the record
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
yd_ind = pd.Index(dt_ind.dayofyear)

# save some info
Info = dict()
Info['run_type'] = Ldir['run_type']
Info['datestring_start'] = dt0.strftime(dsf)
Info['datestring_end'] = dt1.strftime(dsf)

# Load a dataframe with info for rivers to get
gtag = 'cas6_v3'
ri_dir = Ldir['LOo'] / 'pre' / 'river' / gtag
ri_fn = ri_dir / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

# get historical and climatological data files
year0 = 1980
year1 = 2020
# historical and climatological data
Ldir['Hflow_fn'] = ri_dir / 'Data_historical' / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p')
Ldir['Cflow_fn'] = ri_dir / 'Data_historical' / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p')
Ldir['Ctemp_fn'] = ri_dir / 'Data_historical' / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p')

# get the list of rivers and indices for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')

if Ldir['testing']:
    gri_df = gri_df.loc[['columbia'],:]

# associate rivers with ones that have temperature climatology data
ri_df = rivfun.get_tc_rn(ri_df)

# get the flow and temperature data for these days
qt_df_dict = rivfun.get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days)

# get dict S
S_fn = Ldir['grid'] / 'S_COORDINATE_INFO.csv'
S_info_dict = pd.read_csv(S_fn, index_col='ITEMS').to_dict()['VALUES']
S = zrfun.get_S(S_info_dict)

# save the output to NetCDF
out_fn = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc'] / 'rivers.nc'
rivfun.write_to_nc(out_fn, S, gri_df, qt_df_dict, dt_ind)

# add biogeochemical variables
rivfun.add_bio(out_fn, gri_df, yd_ind)

# -------------------------------------------------------

# test for success

if out_fn.is_file():
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
