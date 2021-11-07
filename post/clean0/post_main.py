"""
This is the main program for cleaning out unneeded history files from past forecasts. It
deletes history files 26-73 (or more if the forecast is longer) for the DAY BEFORE the one
it is passed as date_string.

Test on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -job clean0 -test True
Testing just prints what it would do, but does not actually delete the files.

Test on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.05 -job clean0 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -job clean0 > clean.log &

And of course it can be run for a single forecast day or any past range using driver_post.py.
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from lo_tools import Lfun

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)

# the day before date_string
dt1 = dt0 - timedelta(days=1)
ds1 = dt1.strftime(Lfun.ds_fmt)

f_string1 = 'f' + ds1
dir1 = Ldir['roms_out'] / Ldir['gtagex'] / f_string1
if dir1.is_dir():
    print(' - Deleting forecast files for ' + ds1)
    do_delete = True
else:
    print(' - Directory not found for: ' + ds1)
    result_dict['note'] = 'directory not found'
    result = 'fail'
    do_delete = False

keep_name_list = []
for ii in range(1,26):
    iistr = ('0000' + str(ii))[-4:]
    keep_name_list.append('ocean_his_' + iistr + '.nc')
    
if do_delete:
    result = 'success'
    fn_list = [ff for ff in dir1.glob('ocean_his*nc')]
    fn_list.sort()
    jj = 0
    for fn in fn_list:
        if fn.name in keep_name_list:
            pass
        else:
            if Ldir['testing']:
                print('Would delete ' + str(fn))
            else:
                #fn.unlink(missing_ok=True)
                jj += 1
    result_dict['note'] = str(jj) + ' files deleted'

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
