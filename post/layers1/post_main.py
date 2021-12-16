"""
This is the main program for making a box extraction of the daily forecast and sending
the results to the APL server.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -r backfill -job layers1 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's date string] -r forecast -job layers1

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# This will be used for the file name that goes to the public server.
share_name = 'layers'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
out_fn = out_dir / (share_name + '.nc')
temp_dir = out_dir / 'tempfiles'
Lfun.make_dir(temp_dir, clean=True)

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
    fn_list = fn_list[::4]

# do the extractions
N = len(fn_list)
proc_list = []
tt0 = time()
print('Working on ' + Ldir['job'] + ' (' + str(N) + ' times)')
for ii in range(N):
    in_fn = fn_list[ii]
    sys.stdout.flush()
    count_str = ('000000' + str(ii))[-6:]
    temp_out_fn = temp_dir / ('layers_' + count_str + '.nc')
    this_dir = str(Ldir['LO'] / 'post' / Ldir['job']) + '/'
    cmd_list = ['python', this_dir + 'make_layers.py', '-in_fn', str(in_fn),
        '-out_fn', str(temp_out_fn), '-test', str(Ldir['testing'])]
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
    ii += 1
    
# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','layers'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O', str(out_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if True:
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())
print('Time for full layers extraction = %0.2f sec' % (time()- tt0))

if not Ldir['testing']:
    Lfun.make_dir(temp_dir, clean=True)

print('\nPath to file:\n%s' % (str(out_fn)))

# copy the file to the server
post_argfun.copy_to_server(Ldir, out_fn)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
