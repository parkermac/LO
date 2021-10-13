"""
This is the main program for making a LAYERS subset of the daily output.

It creates a single NetCDF file containing fields on selected levels
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job layers0


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

print(' - Creating layers file for ' + Ldir['date_string'])

# this it the name of the file we will copy the output to
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
out_fn = out_dir / 'ocean_layers.nc'

# do the work...

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
