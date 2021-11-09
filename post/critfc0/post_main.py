"""
This is the main program for making the CRITFC output file,
using code from Charles Seaton.

For testing on my mac run in ipython as
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job critfc0 -test True

With -test True it limits the number of files to 3, and prints more to the stdout

Run on apogee
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's date string] -job critfc0 > critfc0.log &

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
import shutil

out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']

# critfc code inputs
rundate = Ldir['date_string'].replace('.','-')
depthfile = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string']) / 'ocean_his_0001.nc'
hgrid = Ldir['data'] / 'critfc' / 'hgrid.ll'
vgrid = Ldir['LO'] / 'post' / Ldir['job'] / 'vgrid.in'
basedir = Ldir['roms_out'] / Ldir['gtagex']
outdir = out_dir

this_dir = str(Ldir['LO'] / 'post' / Ldir['job']) + '/'
cmd = ['python', this_dir + 'gen_cmop_nudge.py', str(hgrid), str(vgrid), str(depthfile),
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

# copy the file to the expected place on boiler
if True:#not Ldir['testing']:
    blr_dir = Path('/boildat/parker/LiveOcean_roms/output/cas6_v3_lo8b/f' + Ldir['date_string'])
    Lfun.make_dir(blr_dir)
    
    cmop_date_string = Ldir['date_string'].replace('.','-')
    salt_name = 'cmop_salt_nu.' + cmop_date_string + '.nc'
    temp_name = 'cmop_temp_nu.' + cmop_date_string + '.nc'
    out_salt_fn = out_dir / salt_name
    out_temp_fn = out_dir / temp_name
    blr_salt_fn = blr_dir / salt_name
    blr_temp_fn = blr_dir / temp_name
    blr_salt_fn.unlink(missing_ok=True)
    blr_temp_fn.unlink(missing_ok=True)
    shutil.copyfile(out_salt_fn, blr_salt_fn)
    shutil.copyfile(out_temp_fn, blr_temp_fn)
    print('\nPath to boiler file:\n%s' % (str(blr_salt_fn)))
    print('Path to boiler file:\n%s' % (str(blr_temp_fn)))
    
    # and then write a little text file to alert the user
    done_fn = blr_dir / 'critfc_done.txt'
    done_fn.unlink(missing_ok=True)
    with open(done_fn, 'w') as ffout:
        ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
    print('Path to done file:\n%s' % (str(done_fn)))

# -------------------------------------------------------

# test for success
if proc.returncode == 0:
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
