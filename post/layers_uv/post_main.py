"""
This is the main program for making an extraction of a limited part of the
model domain.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:
run post_main.py -gtx cas7_t0_x4b -ro 0 -d 2017.07.04 -r backfill -job layers_uv -test True

Run today's forecast on apogee:
python post_main.py -gtx cas6_traps2_x2b -ro 0 -d [today's datestring] -r forecast -job layers_uv > test.log &

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# imports
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import shutil
from lo_tools import Lfun
import numpy as np

print(' - Creating extraction file for ' + Ldir['date_string'])

# this is where the output ends up
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
Lfun.make_dir(out_dir) # usually driver_post1.py would make this, so here we make sure
# it exists for testing.

# create time range for extraction
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
if Ldir['testing']:
    fn_list = fn_list[:2]
else:
    pass

Ncenter = 30
verbose = False
def messages(stdout, stderr, mtitle, verbose):
    # utility function for displaying subprocess info
    if verbose:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(Ncenter,'-'))
            print(stdout.decode())
    if len(stderr) > 0:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        # always print errors
        print(' stderr '.center(Ncenter,'-'))
        print(stderr.decode())
    sys.stdout.flush()

# create the S3 bucket for sharing
fstr = 'f'+ds0
cmd_list = ['s3cmd', 'mb', 's3://'+fstr]
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
messages(stdout, stderr, 's3cmd mb', verbose)

# do the extractions
N = len(fn_list)
proc_list = []
print('Working on ' + Ldir['job'] + ' (' + str(N) + ' times)')
for ii in range(N):
    tt0 = time()
    in_fn = fn_list[ii]
    sys.stdout.flush()
    hour_str = ('000000' + str(ii))[-4:] # name by UTC hour for this forecast day.
    out_fn = out_dir / ('layers_hour_' + hour_str + '.nc')
    this_dir = str(Ldir['LO'] / 'post' / Ldir['job']) + '/'
    cmd_list = ['python', this_dir + 'make_layers.py', '-in_fn', str(in_fn),
        '-out_fn', str(out_fn), '-test', 'False']
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
        proc_list = []
    
if Ldir['testing'] == False:
    for ii in range(N):
        # copy the file to the S3 bucket
        in_fn = fn_list[ii]
        hour_str = ('000000' + str(ii))[-4:] # name by UTC hour for this forecast day.
        out_fn = out_dir / ('layers_hour_' + hour_str + '.nc')
        cmd_list = ['s3cmd', 'put', '--acl-public', str(out_fn), 's3://'+fstr]
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        messages(stdout, stderr, 's3cmd put', verbose)
else:
    pass

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
