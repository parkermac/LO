"""
This is the main program for making the daymovies and pushing them to homer.

Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job daymovie0 -test True

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

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
