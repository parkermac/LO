"""
This is a driver for doing multiple mooring extractions.  It reads in
a dict from LO_user/extract/moor/job_lists.py and uses this to run extract_moor.py as a
series of subprocesses.

2021.11.24 Now it should move the output to a folder named after the job.

Run from the command line like:
python multi_mooring_driver.py -gtx cas6_v3_lo8b -test True > mmd.log &

The same job would be run with flags as:
python multi_mooring_driver.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt hourly -job mickett_2 -get_all True > mmd.log &
NOTE: naming the log as "*.log" means that is it automatically ignored by git (as specifed in LO/.gitignore)
"""

# imports
from lo_tools import Lfun

import sys
import argparse
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
from time import time
import shutil

pid = os.getpid()
print(' multi_mooring_driver '.center(60,'='))
print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly or daily
# select job name
parser.add_argument('-job', type=str) # job name 
# select categories of variables to extract (defined in extract_moor.py)
parser.add_argument('-get_tsa', type=Lfun.boolean_string, default=False)
parser.add_argument('-get_vel', type=Lfun.boolean_string, default=False)
parser.add_argument('-get_bio', type=Lfun.boolean_string, default=False)
parser.add_argument('-get_surfbot', type=Lfun.boolean_string, default=False)
parser.add_argument('-get_pressure', type=Lfun.boolean_string, default=False)
# OR select all of them
parser.add_argument('-get_all', type=Lfun.boolean_string, default=False)
# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
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
# testing
if Ldir['testing']:
    Ldir['roms_out_num'] = 0
    Ldir['ds0'] = '2019.07.04'
    Ldir['ds1'] = '2019.07.06'
    Ldir['list_type'] = 'daily'
    Ldir['job'] = 'scoot'
    Ldir['get_all'] = True
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
# set variable list flags
if Ldir['get_all']:
    Ldir['get_tsa'] = True
    Ldir['get_vel'] = True
    Ldir['get_bio'] = True
    Ldir['get_surfbot'] = True

# get the job_lists module, looking first in LO_user
pth = Ldir['LO'] / 'extract' / 'moor'
upth = Ldir['LOu'] / 'extract' / 'moor'
if (upth / 'job_lists.py').is_file():
    print('Importing job_lists from LO_user')
    job_lists = Lfun.module_from_file('job_lists', upth / 'job_lists.py')
else:
    print('Importing job_lists from LO')
    job_lists = Lfun.module_from_file('job_lists', pth / 'job_lists.py')

# Get job dict:
sta_dict = job_lists.get_sta_dict(Ldir['job'])

# if Ldir['testing']:
#     new_sta_dict = dict()
#     for sn in list(sta_dict.keys())[-5:]:#['F_006_NEW', 'E_006_RIV']:
#         new_sta_dict[sn] = sta_dict[sn]
#     sta_dict = new_sta_dict

# make place for log files
log_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor' / 'logs'
Lfun.make_dir(log_dir)

# make place to copy the results of this job
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor'
jout_dir = out_dir / Ldir['job']
Lfun.make_dir(jout_dir)

print('Results will go to %s' % (str(jout_dir)))

ii = 1
jout_fn_list = []
njobs = len(sta_dict.keys())
for sn in sta_dict.keys():
    tt0 = time()
    print('Working on %s (%d of %d)' % (sn, ii, njobs), end='')
    sys.stdout.flush()
    x = ' ' + str(sta_dict[sn][0])
    y = ' ' + str(sta_dict[sn][1])
    cmd_list = ['python','extract_moor.py',
        '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
        '-0', Ldir['ds0'], '-1', Ldir['ds1'], '-lt', Ldir['list_type'],
        '-sn', sn, '-lon', x, '-lat', y, '-Nproc', str(Ldir['Nproc']),
        '-get_tsa', str(Ldir['get_tsa']), '-get_vel', str(Ldir['get_vel']),
        '-get_bio', str(Ldir['get_bio']), '-get_surfbot', str(Ldir['get_surfbot']),
        '-get_pressure', str(Ldir['get_pressure'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    
    if Ldir['testing']:
        if len(stdout) > 0:
            print(stdout.decode())
        if len(stderr) > 0:
            print(stderr.decode())
        sys.stdout.flush()
    
    # copy or move the results to a folder named for the job
    moor_fn = out_dir / (sn + '_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
    job_moor_fn = jout_dir / (sn + '_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
    try:
        if Ldir['testing']:
            # when testing we keep the original extraction to make it easier to plot
            shutil.copyfile(moor_fn, job_moor_fn)
        else:
            shutil.move(moor_fn, job_moor_fn)
    except FileNotFoundError:
        print(' - error making %s' % (job_moor_fn.name))
    
    # write screen output to logs
    sout_fn = log_dir / (sn + '_screen_output.txt')
    serr_fn = log_dir / (sn + '_subprocess_error.txt')
    sout_fn.unlink(missing_ok=True)
    serr_fn.unlink(missing_ok=True)
    if len(stdout) > 0:
        with open(sout_fn, 'w') as fout:
            fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(serr_fn, 'w') as ffout:
            ffout.write(stderr.decode())
    print(': completed in %d sec' % (time()-tt0))
    sys.stdout.flush()
    ii += 1
print('DONE')