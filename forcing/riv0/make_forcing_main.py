"""
This is the main program for making the RIV forcing file.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -test True -f riv0

"""

from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import forcing_functions as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************
import os
from datetime import timedelta
import pandas as pd

import zrfun
import river_class
import river_functions as rivfun

# set up the time index for the record
dsf = Ldir['ds_fmt']
# set first and last times to be at noon
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
days = (dt0, dt1)
day_list = []
this_day = dt0
while this_day <= dt1:
    day_list.append(this_day)
    this_day += timedelta(days=1)
# pandas Index objects
dt_ind = pd.Index(day_list)
yd_ind = pd.Index(dt_ind.dayofyear)
# save some info
Info = dict()
Info['run_type'] = Ldir['run_type']
Info['datestring_start'] = dt0.strftime(dsf)
Info['datestring_end'] = dt1.strftime(dsf)

#%% Load a dataframe with info for rivers to get
ri_fn = str(Ldir['grid']) + '/river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

#%% associate rivers with ones that have temperature climatology data
df = rivfun.get_tc_rn(df)

# get the flow and temperature data for these days
qt_df_dict = rivfun.get_qt(df, dt_ind, yd_ind, Ldir, dt1, days)

#%% get dict S
S_fn = str(Ldir['grid']) + '/S_COORDINATE_INFO.csv'
S_info_dict = pd.read_csv(S_fn, index_col='ITEMS').to_dict()['VALUES']

S = zrfun.get_S(S_info_dict)

#%% save the output to NetCDF
out_fn = str(Ldir['LOo'] / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']) + '/rivers.nc'
rivfun.write_to_nc(out_fn, S, df, qt_df_dict, dt_ind)

# add biogeochemical variables
rivfun.add_bio(out_fn, df, yd_ind)

# -------------------------------------------------------

# test for success

if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
