"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -test True -f tide0

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import subprocess

cmd = Ldir['which_matlab']

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + date_string) / Ldir['frc']

# This chdir is required to get matlab to find the function if it is
# run by the driver.  If you are running this program on its own it should
# have no effect because you are already in this_dir.
this_dir = str(Path(__file__).absolute().parent)
os.chdir(this_dir)

func = ("make_forcing_worker(\'" +
    str(Ldir['grid']) + "\',\'" +
    str(Ldir['data']) + '/tide' + "\',\'" +
    str(out_dir) + "\',\'" +
    date_string + "\')")
    
if Ldir['testing']:
    print(' function call elements sent to matlab '.center(60,'-'))
    for item in func.split(','):
        print(item)
    
cmd_list = [Ldir['which_matlab'], '-nodisplay', '-r', func, '&']
proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()

# This will end up in Info/screen_output.txt if the code is run by the driver.
print('\n' + ' sdtout '.center(60,'-'))
print(stdout.decode())
print('\n' + ' stderr '.center(60,'-'))
print(stderr.decode())

# -------------------------------------------------------

# test for success
if (len(stderr)==0) and (out_dir / 'tides.nc').is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
