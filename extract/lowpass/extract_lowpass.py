"""
Driver to make tidally averaged files.

This runs over a user-specified range of days for a given gtagex. For each day the
input is 71 hourly files centered on Noon UTC of that day (hence using files from the
day before and the day after). The output is a single NetCDF file called lowpassed.nc
that is placed in the same LO_roms/[gtagex]/[f_string] folder that the central history
files came from. The output file format is very similar to a history file and should
be ready for plotting and extraction just like history files.

This should work even on day 2 of a new run (where we would be missing the 0025 file from
the day before day 1) because we trim the first and last entries from fm_list in lp_worker.py.

To make the job faster this driver splits the work up into Nproc simultaneous subprocesses
which are instances of lp_worker. Each of these works through a subset of the files, multiplying
selected fields by a filter weight and adding them up.

The lp_worker code does its best to automate the selection of fields to low-pass, including
atm fields if they are there and bio fields for old and new code versions if they are there.

Test on mac:
run extract_lowpass -gtx cas7_t0_x4b -0 2017.07.04 -1 2017.07.04 -Nproc 4 -test True

Performance:
mac
-Nproc 4 = 2.5 min per day: BEST CHOICE
-Nproc 10 bogs down my 11-core mac
apogee
-Nproc 10 = 2.6 (check) min per day: BEST CHOICE (16 hours per year)

"""

import sys
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time

from lo_tools import Lfun, zfun, zrfun
from lo_tools import extract_argfun as exfun

Ldir = exfun.intro() # this handles the argument passing

# prepare to loop over all days
ds0 = Ldir['ds0']
ds1 = Ldir['ds1']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

# loop over all days
verbose = 'True'
dtlp = dt0
# dtlp is the day at the middle of the lowpass
while dtlp <= dt1:
    dslp = dtlp.strftime(Lfun.ds_fmt)
    print('\n'+dslp)
    # temporary file output location
    temp_out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'lowpass' / ('temp_' + dslp)
    Lfun.make_dir(temp_out_dir, clean=True)
    # final file output location
    out_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + dslp)

    # Start of chunks loop for this day.
    tt0 = time()
    cca = np.arange(71)
    # Always split the job into Nproc parts. It won't get any faster
    # by using more chunks.
    ccas = np.array_split(cca, Ldir['Nproc'])
    # if Ldir['testing']:
    #     ccas = ccas[:2]
    fnum = 0
    proc_list = []
    N = len(ccas)
    for this_cca in ccas:
        cc0 = this_cca[0]
        cc1 = this_cca[-1]
        # print(' - Working on index range %d:%d' % (cc0, cc1))
        # sys.stdout.flush()
        cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'lowpass' / 'lp_worker.py'),
                    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
                    '-dslp', dslp, '-ii0', str(cc0), '-ii1', str(cc1),
                    '-fnum', str(fnum), '-test', str(Ldir['testing']), '-v', verbose]
        verbose = 'False' # turn off verbose after the first call
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        proc_list.append(proc)
        # Nproc controls how many ncks subprocesses we allow to stack up
        # before we require them all to finish.
        if ((np.mod(fnum,Ldir['Nproc']) == 0) and (fnum > 0)) or (fnum == N-1):
            for proc in proc_list:
                stdout, stderr = proc.communicate()
                if len(stderr) > 0:
                    print(stderr.decode())
                if len(stdout) > 0:
                    print(stdout.decode())
            # make sure everyone is finished before continuing
            proc_list = []
        fnum += 1
    
    # add up the chunks into a single file
    for fnum in range(N):
        fstr = ('000' + str(fnum))[-2:]
        fn = temp_out_dir / ('lp_temp_' + fstr + '.nc')
        ds = xr.open_dataset(fn)
        if fnum == 0:
            lp_full = ds.copy()
        else:
            lp_full = (lp_full + ds).compute()
    # add a time dimension
    lp_full = lp_full.expand_dims('ocean_time')
    lp_full['ocean_time'] = (('ocean_time'), pd.DatetimeIndex([dtlp + timedelta(days=0.5)]))

    # Add other fields to make this easy to use with plotting and extraction tools
    # in the same way as history files. The xarray copy() method returns everything,
    # including dimensions and attributes.
    his_fn = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + ds0) / 'ocean_his_0002.nc'
    ds = xr.open_dataset(his_fn)
    for vn in ['rho','u','v','psi']:
        # note that we have to be explicit for coords
        lp_full.coords['lon_'+vn] = ds.coords['lon_'+vn].copy()
        lp_full.coords['lat_'+vn] = ds.coords['lat_'+vn].copy()
    for vn in ['rho','u','v','psi']:
        lp_full['mask_'+vn] = ds['mask_'+vn].copy()
    for vn in ['s_rho','s_w']:
        lp_full.coords[vn] = ds.coords[vn].copy()
    for vn in ['f','h','hc','Cs_r','Cs_w','Vtransform','pm','pn']:
        lp_full[vn] = ds[vn].copy()

    # also add attributes
    vn_list = list(lp_full.data_vars)
    for vn in vn_list:
        try:
            long_name = ds[vn].attrs['long_name']
            lp_full[vn].attrs['long_name'] = long_name
        except KeyError:
            pass
        try:
            units = ds[vn].attrs['units']
            lp_full[vn].attrs['units'] = units
        except KeyError:
            pass
    ds.close()
    
    if not Ldir['testing']:
        # save output
        out_fn = out_dir / 'lowpassed.nc'
        out_fn.unlink(missing_ok=True)
        lp_full.to_netcdf(out_fn, unlimited_dims=['ocean_time'])
    else:
        out_fn = temp_out_dir / 'lowpassed.nc'
        out_fn.unlink(missing_ok=True)
        lp_full.to_netcdf(out_fn, unlimited_dims=['ocean_time'])
        ds_test = xr.open_dataset(out_fn)
    lp_full.close()

    if not Ldir['testing']:
        # tidying up
        Lfun.make_dir(temp_out_dir, clean=True)
        temp_out_dir.rmdir()
    
    print(' - Time to make tidal average = %0.1f minutes' % ((time()-tt0)/60))
    sys.stdout.flush()
    
    dtlp = dtlp + timedelta(days=1)
