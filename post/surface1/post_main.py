"""
This is the main program for making a SURFACE subset of the daily output.

It creates a single NetCDF file containing only the surface fields
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job surface1


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

print(' - Creating surface file(s) for ' + Ldir['date_string'])

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir / (Ldir['job'] + '_surf_' + Ldir['date_string'] + '_' + Ldir['date_string'] + '.nc')

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ocean_surface.nc'

tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', Ldir['date_string'], '-1', Ldir['date_string'],
    '-lt', 'hourly', '-job', 'surface1', '-surf', 'True', '-uv_to_rho', 'True']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)


stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# copy the file to the expected "post" place
shutil.copyfile(out_fn0, out_fn)

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
