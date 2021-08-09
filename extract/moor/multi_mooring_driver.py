"""
This is a driver for doing multiple mooring extractions.  It reads in
a file from Ldir['data'] and uses this to run extract_moor.py as a
series of subprocesses

Run from the command line like:
python multi_mooring_driver.py -gtx cas6_v3_lo8b -test True > mmd.log &

The same job would be run with flags as:
python multi_mooring_driver.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt hourly -job mickett_2 -get_all True > mmd.log &
NOTE: naming the log as "*.log" means that is it automatically ignored by git (as specifed in LO/.gitignore)
"""

# imports
from pathlib import Path
import sys
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import argparse
import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
from time import time

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
        print('*** Missing required argument to forcing_argfun.intro(): ' + a)
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
    Ldir['roms_out_num'] = 2
    Ldir['ds0'] = '2019.07.04'
    Ldir['ds1'] = '2019.07.06'
    Ldir['list_type'] = 'hourly'
    Ldir['job'] = 'mickett_2'
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

job_dir = Ldir['data'] / 'moor'
if str(job_dir) not in sys.path:
    sys.path.append(str(job_dir))
import job_lists

# Get job dict:
sta_dict = job_lists.get_sta_dict(Ldir['job'])

# make place for log files
log_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor' / 'logs'
Lfun.make_dir(log_dir)

ii = 1
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
        '-sn', sn, '-lon', x, '-lat', y,
        '-get_tsa', str(Ldir['get_tsa']), '-get_vel', str(Ldir['get_vel']),
        '-get_bio', str(Ldir['get_bio']), '-get_surfbot', str(Ldir['get_surfbot'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    
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