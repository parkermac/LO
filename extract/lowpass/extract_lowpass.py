"""
Driver to make tidally averaged files.

It runs over a user-specified range of days for a given gtagex. For each day the
input is 71 hourly files centered on Noon UTC of that day (so using files from the
day before and the day after). The output is a 

Test on mac:
run extract_lowpass -gtx cas7_t0_x4b -0 2017.07.04 -1 2017.07.04 -Nproc 4 -test True

Performance:
mac
-Nproc 4 = 1.7 min per day: BEST CHOICE
-Nproc 10 bogs down my 11-core mac, 2.6 minutes
perigee
-Nproc 20 = 2.5 min per day
-Nproc 10 = 2.0 min per day: BEST CHOICE

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

# get the S dict for future use
S_fn = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + ds0) / 'ocean_his_0002.nc'
S = zrfun.get_basic_info(S_fn, only_S=True)
S_ds = xr.open_dataset(S_fn)
h = S_ds.h
lon_psi = S_ds.lon_psi
lat_psi = S_ds.lat_psi
S_ds.close()

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
    if Ldir['testing']:
        ccas = ccas[:2]
    fnum = 0
    proc_list = []
    N = len(ccas)
    for this_cca in ccas:
        cc0 = this_cca[0]
        cc1 = this_cca[-1]
        # print(' - Working on index range %d:%d' % (cc0, cc1))
        # sys.stdout.flush()
        cmd_list = ['python', 'lp_worker.py',
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
    # add z fields
    # NOTE: this only works if you have h, zeta, and salt as saved fields
    lp_full['h'] = h
    NT, N, NR, NC = lp_full.salt.shape
    lp_full.update({'z_rho':(('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N, NR, NC)))})
    lp_full.update({'z_w':(('ocean_time', 's_w', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N+1, NR, NC)))})
    lp_full.z_rho.attrs = {'units':'m', 'long_name': 'vertical position on s_rho grid, positive up'}
    lp_full.z_w.attrs = {'units':'m', 'long_name': 'vertical position on s_w grid, positive up'}
    hh = h.values
    zeta = lp_full.zeta[0,:,:].values # why is the first index 0? maybe a singletn dimension
    z_rho, z_w = zrfun.get_z(hh, zeta, S)
    lp_full['z_rho'][0,:,:,:] = z_rho
    lp_full['z_w'][0,:,:,:] = z_w
    lp_full.coords['lon_psi'] = (('eta_psi','xi_psi'), lon_psi.values)
    lp_full.coords['lat_psi'] = (('eta_psi','xi_psi'), lat_psi.values)
    out_fn = out_dir / 'lowpassed.nc'
    out_fn.unlink(missing_ok=True)
    lp_full.to_netcdf(out_fn)
    lp_full.close()
    
    # tidying up
    Lfun.make_dir(temp_out_dir, clean=True)
    temp_out_dir.rmdir()
    
    print(' - Time to make tidal average = %0.1f minutes' % ((time()-tt0)/60))
    
    dtlp = dtlp + timedelta(days=1)
