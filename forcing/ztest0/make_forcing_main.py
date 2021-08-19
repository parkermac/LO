"""
This is the main program for making the ZTEST forcing file.  This is intended
to be a template to copy for making all the other make_forcing_main.py
in forcing.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -f ztest0 -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

print(' Ldir seen by make_forcing_main '.center(60,'+'))
for k in Ldir.keys():
    print('%20s : %s' % (k, Ldir[k]))
    
# -------------------------------------------------------

# test for success
if True:
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
