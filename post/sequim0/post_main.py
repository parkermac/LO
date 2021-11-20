"""
This is the main program for making a box extraction of the daily forecast for
Zhaoqing Yang at PNNL.  The region is Sequim Bay


Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job sequim0 -test True

Run for real on apogee:

python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's date string] -job sequim0 > sequim0.log &

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
from lo_tools import Lfun

print(' - Creating extraction file for ' + Ldir['date_string'])

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (Ldir['job'] + '_' + Ldir['date_string'] + '_' + Ldir['date_string'] + '.nc')

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'sequim0.nc'

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', Ldir['date_string'], '-1', Ldir['date_string'],
    '-lt', 'allhours', '-job', 'sequim0']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# copy the file to the expected place on boiler
if not Ldir['testing']:
    blr_dir = Path('/boildat/parker/LiveOcean_roms/output/cas6_v3_lo8b/f' + Ldir['date_string'])
    Lfun.make_dir(blr_dir)
    blr_fn = blr_dir / 'sequim0.nc'
    blr_fn.unlink(missing_ok=True)
    shutil.copyfile(out_fn0, blr_fn)
    print('\nPath to boiler file:\n%s' % (str(blr_fn)))
    
    # and then write a little text file to alert the user
    done_fn = blr_dir / 'sequim_done.txt'
    done_fn.unlink(missing_ok=True)
    with open(done_fn, 'w') as ffout:
        ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
    print('Path to done file:\n%s' % (str(done_fn)))

# move the extraction to the expected "post" place
if Ldir['testing']:
    shutil.copyfile(out_fn0, out_fn)
else:
    shutil.move(out_fn0, out_fn)
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
