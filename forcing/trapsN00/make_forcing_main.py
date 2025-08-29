"""
This is the main program for making the RIVER and TRAPS
forcing file for the updated ROMS

By default, point sources and tiny rivers are enabled. 
To turn them off, set lines 41 and 42 to be = False

Test in ipython:
run make_forcing_main.py -g cas7 -r backfill -d 2019.07.10 -tP trapsP01 -f trapsN00
run make_forcing_main.py -g cas7 -r forecast -d [today's date string] -tP trapsP01 -f trapsN00

where tP is the traps climatology folder
and f is the forcing name (current folder)

NOTES:

I. There are three things you have to do before this can run:

(1) Copy LO_data/trapsD01 (e.g. from /dat1/parker on apogee) into your own LO_data.
(2) Run LO/pre/trapsP01/traps_placement.py once for each grid. I did not have to reprocess any of the raw data in trapsD01 to do this, but you might have to at some point.
(3) Run LO/pre/trapsP01/climatology_all_source_types.sh to create climatologies. To do this you will likely have to do:
  chmod u+x ./climatology_all_source_types.sh
and then
  ./climatology_all_source_types.sh.
(4) THEN you can run the test command for make_forcing_main.py.

II. This code has been modified (2025.08.26 Parker) to write forcing to separate day folders when the run_type = forecast.

"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import sys
from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
from pathlib import Path
from importlib import reload
import rivfun
import make_moh20_LOriv_forcing as LOriv
import make_moh20_triv_forcing as triv
import make_moh20_wwtp_forcing as wwtp_moh20
import make_was24_wwtp_forcing as wwtp_was24

reload(LOriv)
reload(triv)
reload(wwtp_moh20)
reload(wwtp_was24)

#################################################################################
#                     Switch to enable/disable TRAPS                            #
#################################################################################

# LOGICAL SWITCH TO ENABLE OR DISABLE TINY RIVERS OR WWTPS
enable_trivs = True
enable_wwtps = True

#################################################################################
#                              Argument parsing                                 #
#################################################################################

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

#################################################################################
#                    Get required data to generate forcing                      #
#################################################################################

# get correct version of traps climatology output
trapsP = Ldir['trapsP']

# Specify which "ctag" to use. This is related to the collections for which we have
# river and point source climatologies.
ctag = 'lo_base'

# get correct version of Ecology data (based on what is saved in LO/pre/trapsP##)
this_dir = Path(__file__).absolute().parent.parent.parent.parent
with open(this_dir / 'LO' / 'pre' / trapsP / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

if Ldir['testing']:
    reload(zrfun)
    reload(rivfun)

# Set up the time index for the record.
# By default we make it long enough to do a daily forecast.
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=Ldir['forecast_days']+1.5)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
yd_ind = pd.Index(dt_ind.dayofyear)

ot_vec = np.array([Lfun.datetime_to_modtime(item) for item in dt_ind])
NT = len(ot_vec)

S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
N = S['N']

grid_fn = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(grid_fn, only_G=True)

#################################################################################
#                   Run helper scripts to generate forcing                      #
#################################################################################

ds_to_merge = []

# generate forcing for pre-existing LO rivers
try:
    LOriv_ds, NRIV = LOriv.make_forcing(N,NT,dt_ind,yd_ind,ot_vec,dt1,days,Ldir,trapsP,trapsD,ctag)
    ds_to_merge.append(LOriv_ds)
except Exception as e:
    print('Error creating LOriv: maybe there are none.')
    print(e)
    NRIV = 0

# generate forcing for tiny rivers
try:
    triv_ds, NTRIV = triv.make_forcing(N,NT,NRIV,dt_ind,yd_ind,ot_vec,Ldir,enable_trivs,trapsP,trapsD,ctag)
    ds_to_merge.append(triv_ds)
except Exception as e:
    print('Error creating triv: maybe there are none.')
    print(e)
    NTRIV = 0

# generate forcing for marine point sources (Mohamedali et al., 2020)
try:
    wwtp_ds, NWWTP_moh = wwtp_moh20.make_forcing(N,NT,NRIV,NTRIV,dt_ind,yd_ind,ot_vec,Ldir,enable_wwtps,trapsP,trapsD,ctag)
    ds_to_merge.append(wwtp_ds)
except Exception as e:
    print('Error creating moh20 wwtp: maybe there are none.')
    print(e)
    NWWTP_moh = 0

# generate forcing for marine point sources (Wasielewski et al., 2024)
try:
    wwtp_ds, NWWTP_was = wwtp_was24.make_forcing(N,NT,NRIV,NTRIV,NWWTP_moh,dt_ind,yd_ind,ot_vec,Ldir,enable_wwtps,trapsP,trapsD,ctag)
    ds_to_merge.append(wwtp_ds)
except Exception as e:
    print('Error creating was24 wwtp: maybe there are none.')
    print(e)
    NWWTP_was = 0

#################################################################################
#                   Combine forcing outputs and save results                    #
#################################################################################

# combine all forcing datasets
if len(ds_to_merge) == 0:
    print('** No river data at all! **')
    sys.exit()
elif len(ds_to_merge) == 1:
    all_ds = ds_to_merge[0].copy()
else:
    all_ds = xr.merge(ds_to_merge)

# If we are doing a forecast just copy the results to subsequent days.
date_string_list = []
if Ldir['run_type'] == 'backfill':
    date_string_list.append(Ldir['date_string'])
elif Ldir['run_type'] == 'forecast':
            ds00 = Ldir['date_string']
            dt00 = datetime.strptime(ds00, Ldir['ds_fmt'])
            for day in range(Ldir['forecast_days']):
                dt0 = dt00 + timedelta(days=day)
                ds0 = datetime.strftime(dt0, format=Ldir['ds_fmt'])
                date_string_list.append(ds0)
out_fn_list = []
for date_string in date_string_list:
    print('\nWorking on ' + date_string)
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + date_string) / Ldir['frc']
    Lfun.make_dir(out_dir)
    out_fn = out_dir / 'rivers.nc'
    out_fn_list.append(out_fn) # used at the end to check results
    out_fn.unlink(missing_ok=True)
    # Save to NetCDF
    all_ds.to_netcdf(out_fn)
all_ds.close()

# test for success
result = 'SUCCESS'
for fn in out_fn_list:
    if not fn.is_file():
        result = 'FAIL'
result_dict['result'] = result

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
