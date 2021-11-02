"""
This is the main program for reorganizing the files from a daily forecast.
It is designed to work only on a daily forecast folder that has history
files 1-73 in it, and will exit if that is not what it finds.

It creates new folders for the subsequent two days and populates them
with copies of the history files:
day 2: 25-49 become 1-25
day3: 49-73 become 1-25

Test on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -test True > split.log &
Testing just prints what it would do, but does not actually copy the files.

Ru for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] > split.log &
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

start_time = datetime.now()

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=1)
dt2 = dt0 + timedelta(days=2)
ds1 = dt1.strftime(Lfun.ds_fmt)
ds2 = dt2.strftime(Lfun.ds_fmt)

print(' - Distributing forecast files for ' + ds0)
f_string0 = 'f' + ds0
f_string1 = 'f' + ds1
f_string2 = 'f' + ds2

dir0 = Ldir['roms_out'] / Ldir['gtagex'] / f_string0
dir1 = Ldir['roms_out'] / Ldir['gtagex'] / f_string1
dir2 = Ldir['roms_out'] / Ldir['gtagex'] / f_string2

fn_list = [ff for ff in dir0.glob('ocean_his*nc')]
fn_list.sort()
NF = len(fn_list)
if NF != 73:
    print('Expecting 73 history files but only %d found.' % (NF))
    sys.exit()
    
result = 'success'
try:
    Lfun.make_dir(dir1, clean=True)
    Lfun.make_dir(dir2, clean=True)
    for ii in range(1,26):
        h0 = ('0000' + str(ii+24))[-4:]
        h1 = ('0000' + str(ii))[-4:]
        fn0 = dir0 / ('ocean_his_' + h0 + '.nc')
        fn1 = dir1 / ('ocean_his_' + h1 + '.nc')
        if Ldir['testing']:
            print('\n'+str(fn0))
            print(str(fn1))
        else:
            shutil.copyfile(fn0, fn1)
    for ii in range(1,26):
        h0 = ('0000' + str(ii+48))[-4:]
        h2 = ('0000' + str(ii))[-4:]
        fn0 = dir0 / ('ocean_his_' + h0 + '.nc')
        fn2 = dir2 / ('ocean_his_' + h2 + '.nc')
        if Ldir['testing']:
            print('\n'+str(fn0))
            print(str(fn2))
        else:
            shutil.copyfile(fn0, fn2)
except Exception as e:
    result = 'fail'
    print('Error while trying to scopy files:')
    print(e)

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
