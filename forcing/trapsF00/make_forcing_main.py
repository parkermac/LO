"""
This is the main program for making the RIVER and TRAPS
forcing file for the updated ROMS

By default, point sources and tiny rivers are enabled. 
To turn them off, set lines 41 and 42 to be = False

Test on pc in ipython:
run make_forcing_main.py -g cas7 -r backfill -d 2019.07.04 -tP trapsP## -f trapsF##

where tP is the traps climatology folder
and f is the forcing name (current folder)
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
from pathlib import Path
from importlib import reload
import rivfun
import make_LOriv_forcing as LOriv
import make_triv_forcing as triv
import make_wwtp_forcing as wwtp

reload(LOriv)
reload(triv)
reload(wwtp)

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

# get correct version of Ecology data (based on what is saved in LO_output/pre/trapsP##)
this_dir = Path(__file__).absolute().parent.parent.parent
with open(this_dir / 'pre' / trapsP / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

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

# generate forcing for pre-existing LO rivers
LOriv_ds, NRIV = LOriv.make_forcing(N,NT,dt_ind,yd_ind,ot_vec,dt1,days,Ldir,trapsP,trapsD)

# generate forcing for tiny rivers
triv_ds, NTRIV = triv.make_forcing(N,NT,NRIV,dt_ind,yd_ind,ot_vec,Ldir,enable_trivs,trapsP, trapsD)

# generate forcing for marine point sources
wwtp_ds, NWWTP = wwtp.make_forcing(N,NT,NRIV,NTRIV,dt_ind,yd_ind,ot_vec,Ldir,enable_wwtps,trapsP,trapsD)

#################################################################################
#                   Combine forcing outputs and save results                    #
#################################################################################

# combine all forcing datasets
all_ds = xr.merge([LOriv_ds,triv_ds, wwtp_ds])

# Save to NetCDF
all_ds.to_netcdf(out_fn)
all_ds.close()

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
