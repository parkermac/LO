"""
Custom code to copy files to kopah.

Initially created for a request from Jihun Jung at OSU.

I ran this on klone, using:
 conda activate loenv
 python3 custom_copy_to_kopah.py > custom.log &
and it was fine to logout.

Takes about 7 sec per average file, so under an hour per year.
"""

import pandas as pd
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import sys

# step through a sequence of files, copying each from an s3 bucket to a new s3 bucket
# with public access, using s5cmd.

testing = False

if testing:
    tvec = pd.date_range('2024-01-01','2024-01-02')
else:
    tvec = pd.date_range('2024-01-01','2024-12-31')

for t in tvec:
    tstr = t.strftime('%Y.%m.%d')
    in_fn = 's3://liveocean-pmacc/LO_roms/cas7_t1_x11ab/f' + tstr + '/ocean_avg_0001.nc'
    out_dir = 's3://pm-share/jung2024/f' + tstr +'/'

    tt0 = time()
    cmd_list = ['s5cmd','cp','--acl','public-read',in_fn,out_dir]
    if testing:
        print(cmd_list)
    else:
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        # print('stdout:')
        # print(stdout.decode())
        # print('stderr:')
        # print(stderr.decode())
    print('Time to copy: %0.1f sec' % (time()-tt0))
    print("URL = " + out_dir.replace('s3://','https://s3.kopah.uw.edu/') + 'ocean_avg_0001.nc')
    print('')
    sys.stdout.flush()