"""
This is the main program for archiving the daily forecast to another place.
It copies all three days of the forecast, overwriting existing days 1 and 2 from the
previous forecast.

I created this as a standalone program because it is a very niche task.  It is designed as
a post job that runs with the forecast.  I would avoid using it for backfill jobs because it
has the potential to overwrite things you wanted to keep.

Test on mac or apogee:
run post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -r forecast -job archive0 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -r forecast -job archive0 < /dev/null > post.log &
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from subprocess import Popen as PO
from subprocess import PIPE as PI

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Ldir['ds_fmt'])
if Ldir['run_type'] == 'backfill':
    dt1 = dt0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)

result = 'success'
this_dt = dt0
while this_dt <= dt1:
    this_ds = this_dt.strftime(Ldir['ds_fmt'])
    print(' - archiving files for ' + this_ds)
    f_string = 'f' + this_ds
    in_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    out_dir = Path('/pgdat1') / 'parker' / 'LO_roms' / 'cas6_v0_live'
    if Ldir['testing'] == True:
        # For testing we just print the directories
        print('   in_dir: %s' % (str(in_dir)))
        print('   out_dir: %s' % (str(out_dir)))
    else:
        # Copy in_dir to out_dir (overwrite if existing)
        cmd_list = ['scp','-r', str(in_dir), str(out_dir)]
        proc = PO(cmd_list, stdout=PI, stderr=PI)
        stdout, stderr = proc.communicate()
        Ncenter = 30
        if len(stdout) > 0:
            print(' sdtout '.center(Ncenter,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            result = 'fail'
            print(' stderr '.center(Ncenter,'-'))
            print(stderr.decode())
        sys.stdout.flush()
    this_dt += timedelta(days=1)

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
