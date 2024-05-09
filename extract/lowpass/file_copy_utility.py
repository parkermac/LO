"""
Custom code to help Parker copy lowpassed files to his laptop.

** WARNING: I think this code it not ready for general use. It fails too
often even when I run it with an ethernet cable in my UW office.
Running from the bach command line does not help.

Example call, from my mac:
run file_copy_utility -gtx cas7_t0_x4b -0 2013.01.01 -1 2013.01.03

Issue: I have to manually enter my password for each file transfer.
This was easily solved using the usual steps (on my mac):
ssh-keygen
ssh-copy-id parker@apogee.ocean.washington.edu
as described in the readme for LO_roms_user.

"""

import sys

from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from pathlib import Path
from lo_tools import Lfun
from lo_tools import extract_argfun as exfun
from time import time

Ldir = exfun.intro() # this handles the argument passing

# prepare to loop over all days
ds0 = Ldir['ds0']
ds1 = Ldir['ds1']
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

# loop over all days
dt = dt0
# dt is the day at the middle of the lowpass
while dt <= dt1:

    dstr = datetime.strftime(dt, Lfun.ds_fmt)
    print(dstr)
    sys.stdout.flush()

    # Input file
    in_fn = Path('parker@apogee.ocean.washington.edu:/dat1/parker/LO_roms') \
        / Ldir['gtagex'] / ('f' + dstr) / 'lowpassed.nc' 

    # Output file (copy of input)
    out_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + dstr)
    Lfun.make_dir(out_dir)

    tt0 = time()
    cmd_list = ['scp', str(in_fn), str(out_dir)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stderr) > 0:
        print(stderr.decode())
    if len(stdout) > 0:
        print(stdout.decode())
    # Perhaps this would be the place for error checking, but maybe
    # instead we need to use some command other than "scp".
    print(' time to copy = %0.1f seconds' % (time()-tt0))
    sys.stdout.flush()

    dt = dt + timedelta(days=1)
