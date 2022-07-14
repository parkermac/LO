"""
This is the main program for making a box extraction of the daily forecast and sending
the results to the APL server.

This is for Ramsey Harcourt and the TFO-PISCES sound speed project.

Testing on mac:
run post_main.py -gtx cas6_v0_live -ro 0 -d 2019.07.04 -r backfill -job harcourt1 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's date string] -r forecast -job harcourt1

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# +++ These are the main things to edit when changing jobs. +++

# Set the job name for LO/extract/box/extract_box.py, which refers to a choice
# in LO_user/extract/box/job_definitions.py.
box_job = 'harcourt'
# Note that box_job can be different from Ldir['job'] from the command line
# arguments.  Ldir['job'] controls the naming of where we archive the extraction:
# LO_output/post/[gtagex]/[fstring]/Ldir['job']

# This will be used for the file name that goes to the public server.
share_name = 'harcourt'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# imports
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import shutil
from lo_tools import Lfun

print(' - Creating extraction file for ' + Ldir['date_string'])

# create time range for box extraction
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (box_job + '_' + ds0 + '_' + ds1 + '.nc')

# this it the name of the file we will copy/move the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / (share_name + '.nc')

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', ds0, '-1', ds1, '-lt', 'hourly', '-job', box_job]
# note that we don't pass the testing flag, because that would override the
# choices above
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# copy/move the extraction to the expected "post" place
if Ldir['testing']:
    # when testing we keep the original extraction to make it easier to plot
    shutil.copyfile(out_fn0, out_fn)
else:
    shutil.move(out_fn0, out_fn)
print('\nPath to file:\n%s' % (str(out_fn)))

# copy the file to the server
post_argfun.copy_to_server(Ldir, out_fn)

# Extra task: split the file up into hour chunks to make it easier to download
# from the ship.
tt0 = time()
import xarray as xr
ds = xr.open_dataset(out_fn)
NT = len(ds.ocean_time.values)
ds.close()
for hh in range(NT):
    nstr = ('00'+str(hh))[-2:]
    h_out_fn = str(out_dir / 'Data')+'/hour_'+nstr+'.nc'
    cmd_list = ['ncks','-d','ocean_time,%d,%d' % (hh,hh),str(out_fn),
        '-O',h_out_fn]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stderr) > 0:
        print(stderr.decode())
    # copy the file to the server
    post_argfun.copy_to_server(Ldir, h_out_fn, subdir = 'harcourt_hourly')
    
print('Time to split into hours = %0.2f sec' % (time()-tt0))

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
