"""
This is the main program for making a LAYERS subset of the daily output.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job layers0


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

start_time = datetime.now()

print(' - Creating layers file for ' + Ldir['date_string'])

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ocean_layers.nc'

# do the work...

# >>>>>>>>>>>>>>>>>> OLD CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<

# imports
from datetime import datetime
import netCDF4 as nc
import zrfun
import zfun

import layer_fun

import numpy as np
sys.path.append(os.path.abspath('../../plotting/'))
import pfun
from time import time, sleep
from PyCO2SYS import CO2SYS
import seawater as sw
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
import subprocess

start_time = datetime.now()

print('*** Creating layers file for ' + Ldir['date_string'] + ' ***')
f_string = 'f' + Ldir['date_string']

# Create out_dir
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
out_dir = in_dir # output goes to same place as input

# ======== Create an output file for SCOOT, NANOOS and etc. =============================
# performance with subprocess:
# 4 minutes per day (every 4th hour) on my mac

testing = False
verbose = False

# get (subsampled) list of history files to process
fn_list = layer_fun.get_fn_list(in_dir, testing)

# Initialize output file
out_fn = out_dir + 'ocean_layers.nc'
print(' - Writing to: ' + out_fn)
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# NEW parallelize using subprocess

# function defining the subprocess
def run_sub(in_dir, nn, istart, iend, testing=testing, verbose=verbose):
    cmd = ['python', 'make_layers.py', '-in_dir', in_dir, '-nn', str(nn),
            '-istart', str(istart), '-iend', str(iend),
            '-testing',str(testing), '-verbose', str(verbose)]
    proc = subprocess.Popen(cmd)
    return proc

NN = len(fn_list)
if Ldir['lo_env'] == 'pm_mac':
    NP = 4 # max number of processes
else:
    NP = 19 # e.g. on boiler

print('** Creating layers using %d subprocesses for %d files **' % (NP, NN))

# create list of indices of the file list for each job
ii_list = np.array_split(np.arange(NN), NP)

# run a number of jobs
start = time()
procs = []
for pp in range(NP):
    ii = ii_list[pp]
    if len(ii) >= 1:
        istart = ii[0]
        iend = ii[-1]
        print('process %d: istart = %d, iend = %d' % (pp, istart, iend))
        sleep(1) # cludge: needed so that calls to Lstart() don't collide while writing lo_info.csv
        proc = run_sub(in_dir, pp, istart, iend, testing=testing, verbose=verbose)
        procs.append(proc)
    else:
        # ii is an empty array
        pass

# the proc.communicate() method will only return after a job is done
# so this loop effectively checks on all jobs sequentially, and does not
# end until they all are
for proc in procs:
    proc.communicate()

end = time()
print ('Finished in %0.1f sec' % (end - start))

# concatenate to form final file
tlf_list = [in_dir+'temp_layer'+str(ii)+'.nc' for ii in range(NP)]
cmd = ['ncrcat'] + tlf_list + [out_fn]
cproc = subprocess.Popen(cmd)
cproc.communicate()

# clean up
for tlf in tlf_list:
    if os.path.exists(tlf):
      os.remove(tlf)
    else:
      print("The file does not exist")

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
