"""
Code to download WRF forcing files from AtmSci server.

Intended to run on klone, but it is generic.

Intended to be run every morning by cron, just for that day, although
in principle I think you could use it to get whatever days are on the AtmSci server.
"""

import sys, os
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time, sleep
from pathlib import Path

# Add the path to lo_tools by hand so that it we can import Lfun on klone
# without loenv. In general we write code to run on klone using only the
# default python3 installation.
pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

Ldir = Lfun.Lstart()

out_dir0 = Ldir['data'] / 'wrf'
dstr = datetime.now().strftime('%Y%m%d') + '00'
out_dir = out_dir0 / dstr
Lfun.make_dir(out_dir)

print(out_dir)

acct_dict = Lfun.csv_to_dict(Ldir['data'] / 'accounts' / 'wrf_pm_2026.02.15.csv')

cmd_list = ['wget','--no-check-certificate','-e','robots=off','-r','--level=1','--accept','0000','-nd','-nc','-c',
    '--http-user='+acct_dict['username'],
    '--http-passwd='+acct_dict['password'],
    '-P',str(out_dir),
    'https://a.atmos.uw.edu/mm5rt/puget_sound/'+dstr+'/']

tt0 = time()
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stderr) > 0:
    print('Error getting WRF files for %s' % (dstr))
    print(stderr.decode())
print('Time to get WRF files = %0.1f sec' % (time()-tt0))

