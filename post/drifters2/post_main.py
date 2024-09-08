"""
This is the main program for doing a tracker run, converting it
to JSON files, and pushing them to my website for some interactive js.

Testing on mac:

run post_main.py -gtx cas7_t0_x4b -ro 0 -d 2017.07.04 -job drifters2 -test True

Run for real on apogee:

run post_main.py -gtx cas7_t0_x4b -ro 0 -d ['today's datestring] -job drifters2

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from time import time, sleep
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import numpy as np
import json
import xarray as xr
import pandas as pd
from lo_tools import Lfun
import pytz

print('*** Creating drifter files for ' + Ldir['date_string'] + ' ***')
dsr = Ldir['date_string']
if Ldir['testing']:
    dsr0 = dsr
    dtt = 1
else:
    dsr0 = dsr
    dtt = 3

# RUN TRACKER JOBS
tt0 = time()
exp_list = ['PS']
three_flag = 'False'
for exp in exp_list:
    cmd = ['python', str(Ldir['LO']) + '/tracker2/tracker.py',
        '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
        '-d', dsr0, '-dtt', str(dtt),
        '-exp', exp, '-3d', three_flag,
        '-sub_tag', 'forWeb', '-clb', 'True']
    proc = Po(cmd,stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if True:
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
print ('Finished tracker runs in %0.1f sec' % (time() - tt0))
out_fn_dict = {}
out_json_dict0 = {}
out_json_dict1 = {}
for exp in exp_list:
    if three_flag == 'True':
        three_tag = '3d'
    else:
        three_tag = 'surf'
    out_fn_dict[exp] = Ldir['LOo'] / 'tracks2' / Ldir['gtagex'] / (exp + '_' + three_tag + '_forWeb') / ('release_' + dsr0 + '.nc')
    # output location for jsons will be in one folder that is easy to copy, to help with debugging
    json_dir = Ldir['LOo'] / 'tracks2' / Ldir['gtagex'] / 'json_files'
    out_json_dict0[exp] = json_dir / (exp+'_tracks.json')
    out_json_dict1[exp] = json_dir / (exp+'_times.json')
result = 'SUCCESS'
for exp in exp_list:
    if not out_fn_dict[exp].is_file():
        result = 'FAIL'
# END OF TRACKER JOBS

# function to get local time
def get_dt_local(dt, tzl='US/Pacific'):
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

# CONVERT TO JSON AND SCP TO HOMER
for exp in exp_list:

    in_fn = out_fn_dict[exp]
    out_fn0 = out_json_dict0[exp]
    out_fn1 = out_json_dict1[exp]
    Lfun.make_dir(json_dir)
    print('Writing jsons to: ' + str(json_dir))

    ds = xr.open_dataset(in_fn)
    # packed time, particle
    lon = ds['lon'].values
    lat = ds['lat'].values
    NT, NP = lon.shape

    # Make a time vector as a list of strings for the time slider, converted to local time.
    dt0 = ds.ot.values
    dti0 = pd.DatetimeIndex(dt0)
    tt_list = []
    for dtu in dti0:
        dt_local = get_dt_local(dtu)
        tt_list.append(datetime.strftime(dt_local,'%m/%d/%Y - %I%p')+' '+dt_local.tzname())
    
    # Create and save jsons.
    # 1. Position
    xy = []
    for pp in range(NP):
        xy.append({'x': [('%0.3f' % (item)) for item in lon[:,pp]], 'y': [('%0.3f' % (item)) for item in lat[:,pp]]})
    json.dump(xy, open(out_fn0, 'w'))
    # 2. Time
    tt = [{'t': tt_list}]
    json.dump(tt, open(out_fn1, 'w'))
    
    if Ldir['testing']== False:
        
        # send to homer
        cmd2 = ['scp',str(out_json_dict0[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks2/'+exp+'_tracks.json']
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving tracks to homer ' + out_json_dict0[exp].name)
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
            result = 'FAIL'

        cmd2 = ['scp',str(out_json_dict1[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks2/'+exp+'_times.json']
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving times to homer ' + out_json_dict1[exp].name)
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
            result = 'FAIL'
        
    else:
        print('Skipped sending files to homer')
# END CONVERT AND SCP
# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
