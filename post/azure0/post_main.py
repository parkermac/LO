"""
This is the main program for copying the daily forecast to Azure blob storage.
It copies all three days of the forecast, overwriting existing days 1 and 2 from the
previous forecast.

It is designed as a post job that runs with the forecast, but it can also run in backfill.

Test on mac or apogee:
run post_main.py -gtx cas6_v0_live -ro 0 -d 2019.07.04 -r backfill -job azure0 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -r forecast -job azure0 < /dev/null > azure0.log &

Or just put it in the lineup of the driver_post jobs.

Test output on my mac:
====================== copy_to_azure =======================
2022.07.16 11:23:47
 -- Copying: ocean-his-0001.nc to cas6-v0-live-f2019-07-04
 -- URL to access file: https://pm2.blob.core.windows.net/cas6-v0-live-f2019-07-04/ocean-his-0001.nc
 -- Took 183.39 sec
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from time import time
import numpy as np
from lo_tools import Lfun

from subprocess import Popen as Po
from subprocess import PIPE as Pi

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Ldir['ds_fmt'])
if Ldir['run_type'] == 'backfill':
    dt1 = dt0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)

# screen output
print(' copy_to_azure '.center(60,'='))
print(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))

this_dir = Path(__file__).absolute().parent

result = 'success'
this_dt = dt0
while this_dt <= dt1:
    this_ds = this_dt.strftime(Ldir['ds_fmt'])
    f_string = 'f' + this_ds
    fn_list = Lfun.get_fn_list('hourly', Ldir, this_ds, this_ds)

    if Ldir['testing'] == True:
        fn_list = [fn_list[0]]
        
    # Nproc code
    N = len(fn_list)
    proc_list = []
    tt0 = time()
    print('Working on ' + Ldir['job'] + ' (' + str(N) + ' times)')

    for ii in range(N):
        fn = fn_list[ii]
        # AZURE
        # out_name = fn.name.replace('_','-')
        # container_name = Ldir['gtagex'].replace('_','-') + '_' + f_string.replace('.','-')
        out_name = f_string.replace('.','-') + '/' + fn.name.replace('_','-')
        container_name = Ldir['gtagex'].replace('_','-')
        print(' -- Copying: %s to %s' % (out_name, container_name))
        sys.stdout.flush()

        cmd_list = ['python', str(this_dir) + '/forecast_to_azure.py', '-fn', str(fn),
            '-out_name', out_name, '-container_name', container_name]
        if Ldir['testing'] == True:
            print(cmd_list)
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        proc_list.append(proc)
        
        # Nproc controls how many ncks subprocesses we allow to stack up
        # before we require them all to finish.
        # NOTE: we add the (ii > 0) because otherwise it starts by doing a single
        # job, and in this case the jobs are long enough for that to be a significant
        # slowdown.
        if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
            for proc in proc_list:
                stdout, stderr = proc.communicate()
                # make sure everyone is finished before continuing
                if True:
                    if len(stdout) > 0:
                        print('\n'+stdout.decode())
                    if len(stderr) > 0:
                        print('\n'+stderr.decode())
                        result = 'fail'
            proc_list = []
        ii += 1
    print('Time to copy %d files = %0.2f sec' % (N,time()- tt0))
    
    this_dt += timedelta(days=1)

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)

if Ldir['testing'] == True:
    print(result_dict)
