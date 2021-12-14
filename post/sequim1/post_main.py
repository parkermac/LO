"""
This is the main program for making a box extraction of the daily forecast for
Zhaoqing Yang at PNNL.  The region is Sequim Bay


Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -r backfill -job sequim1 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's date string] -r forecast -job sequim1

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

print(' - Creating extraction file for ' + Ldir['date_string'])

# create time range
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)

# name of box job to use (different from Ldir['job'])
box_job = 'sequim0'
share_name = 'sequim' # eventually use this for the public server name

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (box_job + '_' + ds0 + '_' + ds1 + '.nc')

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / (box_job + '.nc')

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', ds0, '-1', ds1,
    '-lt', 'hourly', '-job', box_job]
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# move the extraction to the expected "post" place
if Ldir['testing']:
    shutil.copyfile(out_fn0, out_fn)
else:
    shutil.move(out_fn0, out_fn)
print('\nPath to file:\n%s' % (str(out_fn)))

def messages(stdout, stderr, mtitle, test_flag):
    # utility function for displaying subprocess info
    if test_flag:
        print((' ' + mtitle + ' ').center(60,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
        sys.stdout.flush()

# copy the file to the expected place on the new server
if 'apogee' in Ldir['lo_env']:
    
    share_user = 'parker@liveocean.apl.uw.edu'
    print('Warning: copying extractions to server for sharing only works for parker from apogee.')
    share_dir = '/data/www/liveocean/output/' + 'f' + Ldir['date_string']

    # (i) make the output directory
    cmd_list = ['ssh', share_user, 'mkdir -p ' + share_dir]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    messages(stdout, stderr, 'Make output directory on for sharing', True)
    
    # (ii) copy the extraction to there
    cmd_list = ['scp',str(out_fn), share_user + ':' + share_dir]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    messages(stdout, stderr, 'Copying extraction to ' + share_dir, True)
        
    # (iii) then write a little text file to alert the user
    done_fn = out_fn / (share_name + '_done.txt')
    done_fn.unlink(missing_ok=True)
    with open(done_fn, 'w') as ffout:
        ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
    cmd_list = ['scp',str(done_fn), share_user + ':' + share_dir]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    messages(stdout, stderr, 'Copying done file', True)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
