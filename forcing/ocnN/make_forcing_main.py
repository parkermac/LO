"""
This makes the ocn forcing files for a nested run, interpolating to get
the fields from history files for another run.

Designed to run only as backfill.

Testing:
run make_forcing_main.py -g wgh1 -gtx cas6_traps2_x2b -ro 0 -r backfill -s continuation -d 2017.07.04 -f ocnN -do_bio True -test True

Performance: 6 minutes per day on mac for wgh1 grid with do_bio = True and start_type = new.
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
import numpy as np
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import pickle

from lo_tools import Lfun, zfun, zrfun, Ofun_nc

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# where to find the files to interpolate from
in_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string'])

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# list of history files from the original grid to work from (all Path objects)
h_list = Lfun.get_fn_list('hourly', Ldir, Ldir['date_string'], Ldir['date_string'])

verbose = False
if Ldir['testing']:
    verbose = True
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
         '-his_fn', his_fn, '-start_type', Ldir['start_type'], '-out_fn', str(temp_out_fn),
         '-do_bio', str(Ldir['do_bio'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
        
    if ((np.mod(ii,Nproc) == 0) and (ii > 0)) or (ii == NT-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            # make sure everyone is finished before continuing
            if verbose and (len(stdout) > 0):
                print('\n'+stdout.decode())
            if len(stderr) > 0:
                print('\n'+stderr.decode())
        proc_list = []
    ii += 1
print('Time to run all extractions = %0.1f sec' % (time()-tt0))
sys.stdout.flush()
# +++++++++ end parallel subprocess +++++++++++++++++++++++++++++++++

# Write files to NetCDF.

# Write climatology file making use of zrfun.get_varinfo().
tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
ds = xr.Dataset()
# Initialize Dataset and save variable names to loop over
dd0 = pickle.load(open(temp_out_fn_list[0], 'rb'))
ddk = dd0.keys()
for vn in ddk:
    if vn == 'ocean_time':
        pass
    else:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
        ds[vn] = (dims, np.nan * np.ones((NT,) + dd0[vn].shape))
        ds[vn].attrs['units'] = vinfo['units']
        ds[vn].attrs['long_name'] = vinfo['long_name']
# loop over the temp files
ot_vec = np.nan * np.ones(NT)
ii = 0
for temp_out_fn in temp_out_fn_list:
    dd = pickle.load(open(temp_out_fn, 'rb'))
    for vn in ddk:
        ddv = dd[vn]
        if verbose: # debugging
            print('\n%d %s max = %0.1f' % (ii,vn,np.nanmax(ddv)))
        if vn == 'ocean_time':
            ot_vec[ii] = ddv
        else:
            if ddv.ndim == 2:
                ds[vn][ii, :, :] = ddv
            elif ddv.ndim == 3:
                ds[vn][ii, :, :, :] = ddv
            else:
                print('problem with ndim?')
    ii += 1
# Add the time coordinates
for vn in ddk:
    if vn == 'ocean_time':
        pass
    else:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        tname = vinfo['time_name']
        # time coordinate
        ds[tname] = ((tname,), ot_vec)
        ds[tname].attrs['units'] = Lfun.roms_time_units
# and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
print('- Write clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# Write initial condition file if needed
if Ldir['start_type'] == 'new':
    tt0 = time()
    in_fn = out_dir / 'ocean_clm.nc'
    out_fn = out_dir / 'ocean_ini.nc'
    out_fn.unlink(missing_ok=True)
    Ofun_nc.make_ini_file(in_fn, out_fn)
    print('- Write ini file: %0.2f sec' % (time()-tt0))
    sys.stdout.flush()
else:
    print('- Skipped writing ocean_ini.nc')

# Write boundary file
tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc.make_bry_file(in_fn, out_fn)
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
if Ldir['start_type'] == 'new':
    nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
else:
    nc_list = ['ocean_clm.nc', 'ocean_bry.nc']
if verbose:
    # print info about the files to the screen
    for fn in nc_list:
        print_info(out_dir / fn)
    # check on min depth
    print('')
    ds = xr.open_dataset(out_dir / 'ocean_clm.nc', decode_times=False)
    dsg = xr.open_dataset(grid_fn)
    zz = ds.zeta[0,:,:].values.squeeze()
    hh = dsg.h.values
    print('Minimum Depth = %0.2f m' % (np.nanmin(zz+hh)))
    print('Min zeta = %0.2f m' % (np.nanmin(zz)))
    print('Max zeta = %0.2f m' % (np.nanmax(zz)))
    ds.close()
    dsg.close()
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
