"""
Code to run all the post-processing jobs for a forecast day.

It is designed only to run for a single day.

Testing on mac:

python driver_post.py -gtx cas6_v3_lo8b -r backfill -d 2019.07.04 -ro 2 -test True

"""

import sys
import argparse
from datetime import datetime, timedelta
import subprocess
from time import sleep

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
        print('*** Missing required argument to forcing_argfun.intro(): ' + a)
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

# set day to process
if args.run_type == 'forecast':
    Ldir['date_string'] = datetime.now().strftime(Lfun.ds_fmt)
elif args.run_type == 'backfill':
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
if Ldir['testing']:
    maxcount=3
    sleeptime=1
else:
    maxcount=480
    sleeptime=60
his_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string'])
if Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
else:
    ndays = 1
his_fn_list = []
for his_num in range(1, int(24*ndays) + 2):
    his_string = ('0000' + str(his_num))[-4:]
    his_fn_list.append(his_dir / ('ocean_his_' + his_string + '.nc'))

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
        break
        
    ntries += 1
    if ntries >= maxcount:
        print('Never found all history files.')
        sys.exit()
    else:
        sleep(sleeptime)


# loop over all jobs
job_list = ['surface0']
for job in job_list:
    
    # make clean output directories (often just a place for Info)
    out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / job
    Lfun.make_dir(out_dir, clean=True)
    Lfun.make_dir(out_dir / 'Data')
    Lfun.make_dir(out_dir / 'Info')
    
    j_fn = Ldir['LO'] / 'post' / job / 'post_main.py'
    cmd_list = ['python3', str(j_fn),
                '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
                '-r', args.run_type, '-d', Ldir['date_string'],
                '-job', job, '-test', str(args.testing)]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    with open(out_dir / 'Info' / 'screen_output.txt', 'w') as fout:
        fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(out_dir / 'Info' / 'subprocess_error.txt', 'w') as ffout:
            ffout.write(stderr.decode())
    # this is intended to end up in the log that the cron job makes
    res_fn = out_dir / 'Info' / 'results.txt'
    if res_fn.is_file():
        with open(res_fn, 'r') as fout:
            for line in fout:
                print(line.replace('\n',''))
    else:
        print('ERROR: missing results.txt file')
    
    

