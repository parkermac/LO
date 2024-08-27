"""
Code to run all the post-processing jobs for a forecast day.

It is designed only to run for a single forecast, and expects to find the
history files organized into one-day folders.

Testing on mac:
run driver_post1.py -gtx cas7_t0_x4b -r backfill -d 2017.07.04 -ro 0 -test True

Test on apogee:
python driver_post1.py -gtx cas7_t0_x4b -r forecast -ro 0 -test True < /dev/null > test_post.log &

To run for real on apogee
python driver_post1.py -gtx cas7_t0_x4b -r forecast -ro 0 < /dev/null > post.log &

NOTE: the "< /dev/null" appears to be necessary when running by hand and you stay
logged on because (maybe) the daymovie0 job is somehow expecting standard input,
and when it doesn't get it the job is "Stopped" and you have to use "fg" to start it again.
"""

import sys
import argparse
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time, sleep

from lo_tools import Lfun, zfun

parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
# select time period and frequency
parser.add_argument('-r', '--run_type', type=str)   # backfill or forecast
parser.add_argument('-d', '--date_string', default='', type=str) # e.g. 2019.07.04
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()

gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# set days to process

# first get the start day
if Ldir['run_type'] == 'forecast':
    Ldir['date_string'] = datetime.now().strftime(Lfun.ds_fmt)
elif Ldir['run_type'] == 'backfill':
    if len(Ldir['date_string'])==len(Lfun.ds_fmt)+2:
        pass # assume a valid date was given
    else:
        print('Error: date_string needed for run_type = backfill')
        sys.exit()
else:
    print('Error: Unknown run_type')
    sys.exit()
print((' Post-processing %s for %s' % (Ldir['run_type'], Ldir['date_string'])).center(60,'-'))

# check that all history files are in place
maxcount=1440; sleeptime=60 # will keep looking for 24 hours, e.g. 5 AM to 5 AM the next day
    
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)
his_fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)

if Ldir['testing'] == False:
    all_found = False
    ntries = 0
    while all_found == False:
        for his_fn in his_fn_list:
            if his_fn.is_file():
                all_found = True
                
            elif not his_fn.is_file():
                all_found = False
                break
        if all_found:
            print('All files found. Beginning post-processing.\n')
            sys.stdout.flush()
            sleep(60) # make sure all copying is able to finish
            break
            
        ntries += 1
        if ntries >= maxcount:
            print('Never found all history files.')
            sys.exit()
        else:
            sleep(sleeptime)

tt0 = time()
# loop over all jobs
if Ldir['testing'] == True:
    job_list = ['lowpass0']
else:
    job_list = ['nest_wgh', 'surface1', 'layers1', 'ubc1', 'sequim1',
        'daymovie0', 'drifters0','layers_uv','lowpass0']

for job in job_list:
    
    # make clean output directories (often just a place for Info)
    out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / job
    Lfun.make_dir(out_dir, clean=True)
    Lfun.make_dir(out_dir / 'Data')
    Lfun.make_dir(out_dir / 'Info')
    
    if Ldir['testing']:
        print(job.center(60, '-'))
        print(str(out_dir))
    
    j_fn = Ldir['LO'] / 'post' / job / 'post_main.py'
    cmd_list = ['python3', str(j_fn),
                '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
                '-r', Ldir['run_type'], '-d', Ldir['date_string'],
                '-job', job, '-test', str(Ldir['testing'])]

    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    with open(out_dir / 'Info' / 'screen_output.txt', 'w') as fout:
        fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(out_dir / 'Info' / 'subprocess_error.txt', 'w') as ffout:
            ffout.write(stderr.decode())
            
    # the screen output below is intended to end up in the log that the cron job makes
    res_fn = out_dir / 'Info' / 'results.txt'
    if res_fn.is_file():
        with open(res_fn, 'r') as fout:
            for line in fout:
                print(line.replace('\n',''))
    else:
        print('ERROR: missing results.txt file')
    print('')
    sys.stdout.flush()
    
print('Total time for all post jobs = %0.1f sec' % (time()-tt0))
    
    

