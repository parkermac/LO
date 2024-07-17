"""
This is the main program for making an extraction of a limited part of the
model domain.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:
run post_main.py -gtx cas7_t0_x4b -ro 0 -d 2017.07.04 -r backfill -job layers_uv -test True

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
# Note that in this case we will sand the results to an S3 bucket on kopah.
share_name = 'layers_uv'

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

# create the S3 bucket for sharing
fstr = 'f'+ds0
cmd_list = ['s3cmd', 'mb', 's3://'+fstr]
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
proc.communicate()

# do the extractions
N = len(fn_list)
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
    proc.communicate()
    print('- Hour %s took %0.2f seconds' % (hour_str, time()-tt0))
    sys.stdout.flush()
    
    # copy the file to the S3 bucket
    tt0 = time()
    cmd_list = ['s3cmd', 'put', '--acl-public', str(out_fn), 's3://'+fstr]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc.communicate()
    print('-- copy to S3 bucket %s took %0.2f seconds' % (hour_str, time()-tt0))
    sys.stdout.flush()



# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
