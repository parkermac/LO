"""
This is the main program for making the RIVER and TRAPS
forcing file for the updated ROMS

By default, point sources and tiny rivers are enabled. 
To turn them off, set lines 41 and 42 to be = False

Test in ipython:
run make_forcing_main.py -g cas7 -r backfill -d 2019.07.10 -tP trapsP01 -f trapsN00

where tP is the traps climatology folder
and f is the forcing name (current folder)

2024.10.15 Parker moved the specification of the ctag to a single place, around line 67.
Also added code to catch cases where a grid might not have one or more river types.
These changes were prompted by wanting to use this code for a nested grid.
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

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + date_string) / Ldir['frc']

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

trapsD = 'trapsD01'

if Ldir['testing']:
    reload(zrfun)
    reload(rivfun)

out_fn = out_dir / 'rivers.nc'
out_fn.unlink(missing_ok=True)

# set up the time index for the record
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
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
    print('Error creating wwtp: maybe there are none.')
    print(e)
    NWWTP_moh = 0

# generate forcing for marine point sources (Mohamedali et al., 2020)
try:
    wwtp_ds, NWWTP_was = wwtp_was24.make_forcing(N,NT,NRIV,NTRIV,NWWTP_moh,dt_ind,yd_ind,ot_vec,Ldir,enable_wwtps,trapsD)
    ds_to_merge.append(wwtp_ds)
except Exception as e:
    print('Error creating wwtp: maybe there are none.')
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

# Save to NetCDF
all_ds.to_netcdf(out_fn)
all_ds.close()

# test for success
if out_fn.is_file():
    result_dict['result'] = 'SUCCESS'
else:
    result_dict['result'] = 'FAIL'

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
