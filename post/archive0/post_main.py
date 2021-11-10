"""
This is the main program for archiving the daily forecast to another place.
It just save history files 1-25 for the current day.  It can be run before or after split0.

I created this as a standalone program because it a a very niche task.  It is desigend as
a post job that runs with the forecast, but you could also run it by hand for other days.

NOTE: the target directory must not already exist, or if it does it must be empty.

Test on apogee:
run post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's datestring] -job archive0 -test True
python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's datestring] -job archive0 -test True > archive0.log &
Testing just prints what it would do, but does not actually copy the files.

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's datestring] -job archive0 > archive0.log &
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import shutil
from lo_tools import Lfun

ds0 = Ldir['date_string']

print(' - archiving forecast files for ' + ds0)
f_string = 'f' + ds0

in_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
out_dir = Path('/pgdat1') / 'parker' / 'LiveOcean_roms' / 'output' / 'cas6_v3_lo8b' / f_string

if out_dir.is_dir():
    try:
        out_dir.rmdir()
    except Exception as e:
        print('Error - target directory exists and is not empty!')
        print(e)
        sys.exit()

name_list = []
for ii in range(1,26):
    hh = ('0000' + str(ii))[-4:]
    name_list.append('ocean_his_' + hh + '.nc')

result = 'success'
do_copy = True
for name in name_list:
    fn = in_dir / name
    if fn.is_file():
        pass
    else:
        do_copy = False
        result = 'fail'
        result_dict['note'] = 'Missing file ' + name
        break
    
result = 'success'
if do_copy:
    try:
        Lfun.make_dir(out_dir)
        for name in name_list:
            in_fn = in_dir / name
            out_fn = out_dir / name
            if Ldir['testing']:
                print('\n'+str(in_fn))
                print(str(out_fn))
            else:
                shutil.copyfile(in_fn, out_fn)
    except Exception as e:
        result = 'fail'
        print('Error while trying to copy files:')
        print(e)

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
