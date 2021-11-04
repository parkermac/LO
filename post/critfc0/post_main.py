"""
This is the main program for making the CRITFC output file,
using code from Charles Seaton.

For testing on my mac run in ipython as
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job critfc0 -test True

With -test True it limits the number of files to 3, and prints more to the stdout

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time

out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']

# critfc code inputs
rundate = Ldir['date_string'].replace('.','-')
depthfile = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string']) / 'ocean_his_0001.nc'
hgrid = Ldir['data'] / 'critfc' / 'hgrid.ll'
vgrid = Ldir['LO'] / 'post' / Ldir['job'] / 'vgrid.in'
basedir = Ldir['roms_out'] / Ldir['gtagex']
outdir = out_dir

cmd = ['python', 'gen_cmop_nudge.py', str(hgrid), str(vgrid), str(depthfile),
    str(basedir), str(outdir), rundate, '-test', str(Ldir['testing'])]
proc = Po(cmd, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if True:
    if len(stdout) > 0:
        print(' sdtout '.center(60,'-'))
        print(stdout.decode())
    if len(stderr) > 0:
        print(' stderr '.center(60,'-'))
        print(stderr.decode())

# -------------------------------------------------------

# test for success
if proc.returncode == 0:
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
