"""
This makes the ocn forcing files for a nested run, interpolating to get
the fields from history files for another run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g hc0 -t v0 -r backfill -s continuation -d 2019.07.04 -f ocnN -gtx_nest cas6_traps2_x2b -ro_nest 0 -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
from scipy.spatial import cKDTree
import numpy as np
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import pickle

from lo_tools import Lfun, zfun, zrfun

import Ofun_nc_xarray
# defaults
verbose = True
if Ldir['testing']:
    # verbose = True
    from importlib import reload
    reload(Ofun_nc_xarray)

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# where to find the files to interpolate from
in_dir = Ldir['roms_out_nest'] / Ldir['gtagex_nest'] / ('f' + Ldir['date_string'])

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# list of history files from the original grid to work from (all Path objects)
h_list = sorted(in_dir.glob('ocean_his_*'))
# NEW
h_list = Lfun.get_fn_list('hourly', Ldir, Ldir['date_string'], Ldir['date_string'])
if Ldir['testing']:
    h_list = h_list[:2]

# +++++++++++ parallel subprocess +++++++++++++++++++++++++++++++++++
temp_dir = out_dir / 'Data'
Lfun.make_dir(temp_dir, clean=True)
Nproc = 10
NT = len(h_list)
proc_list = []
tt0 = time()
temp_out_fn_list = []
print('Working on ' + Ldir['frc'] + ' (' + str(NT) + ' times)')
for ii in range(NT):
    grid_fn = str(Ldir['grid'] / 'grid.nc')
    his_fn = str(h_list[ii])
    count_str = ('000000' + str(ii))[-6:]
    temp_out_fn = temp_dir / ('data_dict_' + count_str + '.p')
    temp_out_fn_list.append(temp_out_fn)
    this_dir = str(Ldir['LO'] / 'forcing' / Ldir['frc']) + '/'
    cmd_list = ['python', this_dir + 'get_one_time.py', '-grid_fn', grid_fn,
         '-his_fn', his_fn, '-start_type', Ldir['start_type'], '-out_fn', str(temp_out_fn)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
        
    if ((np.mod(ii,Nproc) == 0) and (ii > 0)) or (ii == NT-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            # make sure everyone is finished before continuing
            if True:
                if len(stdout) > 0:
                    print('\n'+stdout.decode())
                if len(stderr) > 0:
                    print('\n'+stderr.decode())
        proc_list = []
    ii += 1
print('Time to run all extractions = %0.1f sec' % (time()-tt0))
sys.stdout.flush()
# +++++++++ end parallel subprocess +++++++++++++++++++++++++++++++++

# Write files to NetCDF.

tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_clm_file(temp_out_fn_list, NT, out_fn)
print('- Write clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_ini.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_ini_file(in_fn, out_fn)
print('- Write ini file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_bry_file(in_fn, out_fn)
print('- Write bry file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# clean up
Lfun.make_dir(temp_dir, clean=True)
    
def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn, decode_times=False)
    print(ds)
    ds.close()

# check results
nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
if False:
    # print info about the files to the screen
    for fn in nc_list:
        print_info(out_dir / fn)
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
