"""
This is the main program for making a tidally averaged version of the daily forecast.

It just does one day, the first day of the forecast. It could in principle do the 
second day, but why bother as it will just be overwritten.

Test on apogee:
run post_main.py -gtx cas7_t0_x4b -ro 0 -d [today's datestring] -r forecast -job lowpass0
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from subprocess import Popen as PO
from subprocess import PIPE as PI
ds0 = Ldir['date_string']
result = 'SUCCESS'
print(' - creating lowpassed file for ' + ds0)
f_string = 'f' + ds0
# command to emulate
# run extract_lowpass -gtx cas7_t0_x4b -0 2017.07.04 -1 2017.07.04 -Nproc 10
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'lowpass' / 'extract_lowpass.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', ds0, '-1', ds0, '-Nproc', '10']
# Note that we make use of the cast that testing = False by default. Hence when we run this
# using driver_post1.py with testing = True it will still run normally.
proc = PO(cmd_list, stdout=PI, stderr=PI)
stdout, stderr = proc.communicate()
Ncenter = 30
if len(stdout) > 0:
    print(' sdtout '.center(Ncenter,'-'))
    print(stdout.decode())
if len(stderr) > 0:
    result = 'FAIL'
    print(' stderr '.center(Ncenter,'-'))
    print(stderr.decode())
sys.stdout.flush()

# -------------------------------------------------------

# test for success
in_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
out_fn = in_dir / 'lowpassed.nc'
if not out_fn.is_file():
    result = 'FAIL'
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
