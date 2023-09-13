"""
Code to run one or more post-processing jobs for a range of days in the past.

If -test False, and it is parker running it on apogee, it sends the output
to the APL server.

>>> It only runs as backfill. <<<

Testing on mac:
run driver_post_backfill.py -gtx cas6_v0_live -0 2019.07.04 -ro 0 -test True
Performance: takes 355 sec per day on my mac for layers1.

Test on apogee:
python driver_post_backfill.py -gtx cas6_v0_live -ro 1 -0 2021.01.01 -test True < /dev/null > test_post_backfill.log &

Run for real on apogee
python driver_post_backfill.py -gtx cas6_v0_live -ro 1 -0 2021.01.01 -1 2021.01.02 < /dev/null > post_backfill.log &


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
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v0_live
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex','ds0']:
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
Ldir['run_type'] = 'backfill'
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# set time range to process
ds0 = Ldir['ds0']
if len(Ldir['ds1']) == 0:
    ds1 = ds0
else:
    ds1 = Ldir['ds1']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
    
print(('Post backfill for: %s-%s ' % (ds0, ds1)).center(60,'-'))

# specify jobs to run
if Ldir['testing'] == True:
    job_list = ['layers1']
else:
    job_list = ['layers1']

# loop over all days
dt = dt0
tt0 = time()
while dt <= dt1:
    tt00 = time()
    
    Ldir['date_string'] = datetime.strftime(dt, Lfun.ds_fmt)
    
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
        
    dt += timedelta(days=1)
    print(' - Time for post jobs for one day = %0.1f sec' % (time()-tt00))
    
print('Total time for all post jobs = %0.1f sec' % (time()-tt0))
    
    

