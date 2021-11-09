"""
This is the main program for making a LAYERS subset of the daily output.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job layers0 -test True


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
from lo_tools import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import numpy as np

print('Creating layers file for ' + Ldir['date_string'])

out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ocean_layers.nc'
temp_dir = out_dir / 'tempfiles'
Lfun.make_dir(temp_dir, clean=True)

in_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string'])
fn_list = Lfun.get_fn_list('allhours', Ldir, Ldir['date_string'], Ldir['date_string'])
if Ldir['testing']:
    fn_list = fn_list[:2]
else:
    fn_list = fn_list[::4]
    
# do the extractions
N = len(fn_list)
proc_list = []
tt0 = time()
print('Working on ' + Ldir['job'] + ' (' + str(N) + ' times)')
for ii in range(N):
    in_fn = fn_list[ii]
    sys.stdout.flush()
    count_str = ('000000' + str(ii))[-6:]
    temp_out_fn = temp_dir / ('layers_' + count_str + '.nc')
    this_dir = str(Ldir['LO'] / 'post' / Ldir['job']) + '/'
    cmd_list = ['python', this_dir + 'make_layers.py', '-in_fn', str(in_fn),
        '-out_fn', str(temp_out_fn), '-test', str(Ldir['testing'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
        
    # Nproc controls how many ncks subprocesses we allow to stack up
    # before we require them all to finish.
    # NOTE: we add the (ii > 0) because otherwise it starts by doing a single
    # job, and in this case the jobs are long enough for that to be a significant
    # slowdown.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            # make sure everyone is finished before continuing
            if True:
                if len(stdout) > 0:
                    print('\n'+stdout.decode())
                if len(stderr) > 0:
                    print('\n'+stderr.decode())
        proc_list = []
    ii += 1
    
# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','layers'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O', str(out_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if True:
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())
print('Time for full layers extraction = %0.2f sec' % (time()- tt0))

# copy the file to the expected place on boiler
if not Ldir['testing']:
    blr_dir = Path('/boildat/parker/LiveOcean_roms/output/cas6_v3_lo8b/f' + Ldir['date_string'])
    Lfun.make_dir(blr_dir)
    blr_fn = blr_dir / 'ocean_layers.nc'
    blr_fn.unlink(missing_ok=True)
    shutil.copyfile(out_fn, blr_fn)
    print('\nPath to boiler file:\n%s' % (str(blr_fn)))
    
    # and then write a little text file to alert the user
    done_fn = blr_dir / 'layers_done.txt'
    done_fn.unlink(missing_ok=True)
    with open(done_fn, 'w') as ffout:
        ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
    print('Path to done file:\n%s' % (str(done_fn)))

if not Ldir['testing']:
    Lfun.make_dir(temp_dir, clean=True)

print('\nPath to file:\n%s' % (str(out_fn)))

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
