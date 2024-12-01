"""
This is the main program for making creating the forcing for a nested model
daily forecast.

This is specific to the wgh (Willapa Bay - Grays Harbor) grid.

Testing on mac:

run post_main.py -gtx cas6_traps2_x2b -d 2017.07.04 -r backfill -job nest_wgh
- works for ocnN but fails for atm00

run post_main.py -gtx cas6_traps2_x2b -d 2019.07.04 -r backfill -job nest_wgh
- works for atm00 but fails for ocnN

Run for real on apogee:
python post_main.py -gtx cas6_traps2_x2b -d [today's date string] -r forecast -job nest_wgh > nest_test.log &

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

nest_gridname = 'wgh2'

# imports
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import shutil
from lo_tools import Lfun

print(' - Creating ' + nest_gridname + ' forcing for ' + Ldir['date_string'])

# create time range for box extraction
ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
if Ldir['run_type'] == 'backfill':
    ds1 = ds0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)
    ds1 = dt1.strftime(Lfun.ds_fmt)

"""
Example commands to emulate here:
python driver_forcing3.py -g wgh2 -gtx cas7_t0_x4b -ro 0 -do_bio True -r backfill -s perfect -0 2024.12.01 -1 2024.12.03 -f ocnN > ocnN_wgh2.log &
python driver_forcing3.py -g wgh2 -r forecast -f atm00 > atm00_wgh2.log &
python driver_forcing3.py -g wgh2 -r forecast -f riv00 > riv00_wgh2.log &    
"""
# make forcing

tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'driver' / 'driver_forcing3.py'), '-g', nest_gridname,
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-do_bio', 'True', '-r', 'backfill', '-s', 'perfect',
    '-0', ds0, '-1', ds1, '-f', 'ocnN']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
    result_dict['note'] = 'ocnN problem'
print('Elapsed time = %0.2f sec' % (time()-tt0))

for frc in ['atm00', 'riv00']:
    tt0 = time()
    if Ldir['run_type'] == 'forecast':
        cmd_list = ['python', str(Ldir['LO'] / 'driver' / 'driver_forcing3.py'), '-g', nest_gridname,
            '-r', 'forecast', '-f', frc]
    elif Ldir['run_type'] == 'backfill':
        cmd_list = ['python', str(Ldir['LO'] / 'driver' / 'driver_forcing3.py'), '-g', nest_gridname,
            '-r', 'backfill', '-0', ds0, '-f', frc]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    print(stdout.decode())
    if len(stderr) > 0:
        print(stderr.decode())
        if 'note' in result_dict.keys():
            result_dict['note'] += ', ' + frc + ' problem'
        else:
            result_dict['note'] = frc + ' problem'
    print('Elapsed time = %0.2f sec' % (time()-tt0))


# -------------------------------------------------------

# test for success
if 'note' not in result_dict.keys():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
