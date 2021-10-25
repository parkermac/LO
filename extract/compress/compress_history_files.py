"""
Code to compress, in-place, a sequence of ROMS history files.

This would typically only be applied to older model output, such as from cas6_v3_lo8b,
that is uncompressed, NetCDF 3.

Test on mac:

run compress_history_files -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04
python compress_history_files.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04

Performance: took 3 minutes to compress 1 day of cas6 on my laptop.  This is 18 hours per year,
however that was with Nproc = 4.  With Nproc = 10 on perigee it took 8 hours for a year of
history files with no NPZD.  Using Nproc = 100 did not speed things up, but it did make
perigee a bit sluggish on the command line.  Overall I think 5 sec per file might be typical
(meaning that I am running 10 files at once and each takes 50 sec).

"""

import sys
import argparse
from lo_tools import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import numpy as np
import xarray as xr
from time import time
from datetime import datetime, timedelta

print(' compress history files '.center(60,'='))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
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

# we build fn_list one day at a time to make sure we get all files
dt0 = datetime.strptime(Ldir['ds0'], Lfun.ds_fmt)
dt1 = datetime.strptime(Ldir['ds1'], Lfun.ds_fmt)
ds_list = Lfun.date_list_utility(dt0, dt1)
fn_list = []
for this_ds in ds_list:
    fn_list = fn_list + Lfun.get_fn_list('allhours', Ldir, this_ds, this_ds)

# do the compression
N = len(fn_list)
proc_list = []
tt0 = time()
for ii in range(N):
    fn = fn_list[ii]
    cmd_list = ['python', 'compress_a_file.py', str(fn)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    # complete Nproc jobs before filling proc_list again
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            # make sure everyone is finished before continuing
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                print(stdout.decode().replace('\n',''))
            if len(stderr) > 0:
                print(stderr.decode())
            sys.stdout.flush()
        proc_list = []
    ii += 1
print('Total time to compress %d files = %0.1f sec' % (N, time()-tt0))
    
