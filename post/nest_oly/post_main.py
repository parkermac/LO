"""
This is the main program for making creating the forcing for a nested model
daily forecast.

This is specific to the South Sound grid.

Testing on mac:

run post_main.py -gtx cas7_t0_x4b -d 2017.07.04 -r backfill -job nest_oly
- works for ocnN but fails for atm00

run post_main.py -gtx cas7_t0_x4b -d 2019.07.04 -r backfill -job nest_oly
- works for atm00 but fails for ocnN

Run for real on apogee:
python post_main.py -gtx cas7_t0_x4b -d [today's date string] -r forecast -job nest_oly > nest_test.log &

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

nest_gridname = 'oly1'

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

One day from scratch:
python driver_forcing3.py -g oly1 -gtx cas7_t0_x4b -ro 0 -do_bio True -r backfill -s new -0 2024.10.07 -1 2024.10.07 -f ocnN > ocnN_oly1.log &
python driver_forcing3.py -g oly1 -r backfill -0 2024.10.07 -1 2024.10.07 -f atm00 > atm00_oly1.log &
python driver_forcing3.py -g oly1 -r backfill -0 2024.10.07 -1 2024.10.07 -tP trapsP00 -f trapsF00 > trapsF00_oly1.log &

or
Forecast:
python driver_forcing3.py -g oly1 -gtx cas7_t0_x4b -ro 0 -do_bio True -r backfill -s perfect -0 2024.10.07 -1 2024.10.09 -f ocnN > ocnN_oly1.log &
python driver_forcing3.py -g oly1 -r forecast -f atm00 > atm00_oly1.log &
python driver_forcing3.py -g oly1 -r forecast -tP trapsP00 -f trapsF00 > trapsF00_oly1.log &    
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

frc = 'atm00'
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

frc = 'trapsF00'
tt0 = time()
if Ldir['run_type'] == 'forecast':
    cmd_list = ['python', str(Ldir['LO'] / 'driver' / 'driver_forcing3.py'), '-g', nest_gridname,
        '-r', 'forecast', '-f', frc, '-tp', 'trapsP00']
elif Ldir['run_type'] == 'backfill':
    cmd_list = ['python', str(Ldir['LO'] / 'driver' / 'driver_forcing3.py'), '-g', nest_gridname,
        '-r', 'backfill', '-0', ds0, '-f', frc, '-tP', 'trapsP00']
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
    result_dict['result'] = 'SUCCESS'
else:
   result_dict['result'] = 'FAIL'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
