"""
This is the main program for making a SURFACE subset of the daily output.

It creates a single NetCDF file containing only the surface fields
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job surface0 -test True

Run for real on apogee:

python post_main.py -gtx cas6_v0_u0mb -ro 0 -d [today's date string] -job surface0 > surface.log &

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

print(' - Creating surface file for ' + Ldir['date_string'])

# this is the name of the file created by extract/box/extract_box.py
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
out_fn0 = out_dir0 / (Ldir['job'] + '_surf_' + Ldir['date_string'] + '_' + Ldir['date_string'] + '.nc')

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ocean_surface.nc'

# run extract_box.py to do the actual job
tt0 = time()
cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'box' / 'extract_box.py'),
    '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
    '-0', Ldir['date_string'], '-1', Ldir['date_string'],
    '-lt', 'allhours', '-job', 'surface0', '-surf', 'True', '-uv_to_rho', 'True']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
print(stdout.decode())
if len(stderr) > 0:
    print(stderr.decode())
print('Elapsed time = %0.2f sec' % (time()-tt0))

# Send the file to Azure, to be picked up by IOOS EDS Viewer folks
if not Ldir['testing']:
    f_string = 'f' + Ldir['date_string']
    ff_string = f_string.replace('.','') # azure does not like dots in container names
    container_name = ff_string
    cmd_list = ['python', str(Ldir['LO'] / 'misc' / 'copy_to_azure.py'),
        '-fn', str(out_fn0), '-out_name', 'ocean_surface.nc',
        '-container_name', container_name]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    print(stdout.decode())
    if len(stderr) > 0:
        print(stderr.decode())
else:
    print('<< Skipped copying file to Azure >>')

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
