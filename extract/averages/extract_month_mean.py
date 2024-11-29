"""
Driver to monthly mean files.

This makes monthly averages of lowpassed.nc files.

Where the output goes:
LO_roms/[gtagex]/averages/month_mean_[year]_[month].nc

The file format should be much like history files, and easily plotted and maniputated.

Test on mac:
run extract_month_mean -gtx cas7_t0_x4b -0 2020.01.01 -1 2020.03.31 -Nproc 4 -test True
Just outputs info about input and output files.

Run on apogee:
run extract_month_mean -gtx cas7_t0_x4b -0 2020.01.01 -1 2020.01.31

NOTE: The 0 and 1 inputs just provide the year and month of the start and end for averages.
The days are ignored. Using "-0 2017.01.01 -1 2017.01.31" will just produce
month_mean_2017_01.nc

Performance:
mac
-Nproc ## min per day: BEST CHOICE
-Nproc 10 = ## min per day: BEST CHOICE (## hours per year)

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

dti = pd.date_range(dt0, dt1, freq='ME', inclusive='right')
# This is a DatetimeIndex of the last day of each month in the range,
# including the last month.

dt_dict = dict()
# Now make a dict filled with DatetimeIndex objects of days for each month
for dt in dti:
    dt00 = datetime(dt.year,dt.month,1)
    dt11 = dt
    this_dti = pd.date_range(dt00, dt11, freq='D', inclusive='both')
    dt_dict[dt] = this_dti

# check results (Result: looks good)
if False:
    for dt in dti:
        this_dti = dt_dict[dt]
        print('')
        print(this_dti[0].strftime(Lfun.ds_fmt))
        print(this_dti[-1].strftime(Lfun.ds_fmt))

# final file output location
out_dir = Ldir['roms_out'] / Ldir['gtagex'] / 'averages'
Lfun.make_dir(out_dir)

# loop over all months
verbose = 'True'
for dt in dti:

    this_ym = dt.strftime('%Y_%m')
    print('\n'+this_ym)

    # temporary file output location
    temp_out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'averages' / ('temp_' + this_ym)
    if verbose:
        print('temp_out_dir')
        print(str(temp_out_dir))
    Lfun.make_dir(temp_out_dir, clean=True)

    # Start of chunks loop for this day.
    tt0 = time()

    this_dti = dt_dict[dt]
    days_per_month = len(this_dti)
    cca = np.arange(days_per_month)
    # Always split the job into Nproc parts. It won't get any faster
    # by using more chunks.
    ccas = np.array_split(cca, Ldir['Nproc'])
    fnum = 0
    proc_list = []
    N = len(ccas)
    for this_cca in ccas:
        this_ds0 = this_dti[this_cca[0]].strftime(Lfun.ds_fmt)
        this_ds1 = this_dti[this_cca[-1]].strftime(Lfun.ds_fmt)
        if verbose:
            print('  fnum = ' + str(fnum))
            print('   ' + this_ds0)
            print('   ' + this_ds1)
            sys.stdout.flush()

        cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'averages' / 'month_mean_worker.py'),
                    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
                    '-this_ym', this_ym, '-0', this_ds0, '-1', this_ds1,
                    '-days_per_month', str(days_per_month),
                    '-fnum', str(fnum), '-test', str(Ldir['testing']), '-v', verbose]
        if verbose:
            print('   ' + (' ').join(cmd_list))
        
        if not Ldir['testing']:
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
        fn = temp_out_dir / ('mean_temp_' + fstr + '.nc')

        if Ldir['testing']:
            print('  + ' + str(fn))
        else:
            ds = xr.open_dataset(fn, decode_times=False)
            if fnum == 0:
                mean_full = ds.copy()
            else:
                mean_full = (mean_full + ds).compute()
    
    # save output
    out_fn = out_dir / ('monthly_mean_' + this_ym + '.nc')

    if Ldir['testing']:
        print(str(out_fn))
    else:
        out_fn.unlink(missing_ok=True)
        mean_full.to_netcdf(out_fn, unlimited_dims=['ocean_time'])

        if Ldir['testing']:
            ds_test = xr.open_dataset(out_fn)
        mean_full.close()

    if not Ldir['testing']:
        # tidying up
        Lfun.make_dir(temp_out_dir, clean=True)
        temp_out_dir.rmdir()
    
    print(' - Time to make monthly_mean = %0.1f minutes' % ((time()-tt0)/60))
    sys.stdout.flush()

