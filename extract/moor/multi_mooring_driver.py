"""
This is a driver for doing multiple mooring extractions.  It reads in
a file from Ldir['data'] and uses this to run extract_moor.py as a
series of subprocesses

Run from the command like:

python multi_mooring_driver.py -job mickett_2 -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -get_tsa True -get_vel True -get_bio True > dlog &

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun

Ldir = exfun.intro() # this handles the argument passing

job_dir = Ldir['data'] / 'moor'
if str(job_dir) not in sys.path:
    sys.path.append(str(job_dir))
import job_lists

import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi

# Get job dict:
sta_dict = job_lists.get_sta_dict(Ldir['job'])

log_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor' / 'logs'
Lfun.make_dir(log_dir)

for sn in sta_dict.keys():
    x = ' ' + str(sta_dict[sn][0])
    y = ' ' + str(sta_dict[sn][1])
    cmd_list = ['python','extract_moor.py',
        '-g', Ldir['gridname'], '-t', Ldir['tag'], '-x', Ldir['ex_name'], '-ro', str(Ldir['roms_out_num']),
        '-0', Ldir['ds0'], '-1', Ldir['ds1'],
        '-sn', sn, '-lon', x, '-lat', y,
        '-get_tsa', str(Ldir['get_tsa']), '-get_vel', str(Ldir['get_vel']), '-get_bio', str(Ldir['get_bio'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    
    # write screen output to logs
    sout_fn = log_dir / (sn + '_screen_output.txt')
    serr_fn = log_dir / (sn + '_subprocess_error.txt')
    sout_fn.unlink(missing_ok=True)
    serr_fn.unlink(missing_ok=True)
    if len(stdout) > 0:
        with open(sout_fn, 'w') as fout:
            fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(serr_fn, 'w') as ffout:
            ffout.write(stderr.decode())
