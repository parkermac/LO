"""
This is the main program for doing a tracker run, converting it
to JSON files, and pushing them to my website for some interactive js.

Testing on mac:

run post_main.py -gtx wgh2_t0_xn0b -ro 0 -d 2024.08.12 -job drifters_wgh -test True

Run for real on apogee:

run post_main.py -gtx wgh2_t0_xn0b -ro 0 -d ['today's datestring] -job drifters_wgh

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

print('*** Creating drifter files for ' + Ldir['date_string'] + ' ***')
dsr = Ldir['date_string']
if Ldir['testing']:
    dtt = 1
else:
    dtt = 3 # assume this is for a daily forecast, so we track for three days

# RUN TRACKER JOBS - parallelize using subprocess
tt0 = time()
procs = []
exp_list = ['wgh0']
three_flag = 'True'
for exp in exp_list:
        cmd = ['python', str(Ldir['LO']) + '/tracker2/tracker.py',
            '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
            '-d', dsr, '-dtt', str(dtt),
            '-exp', exp, '-3d', three_flag,
            '-sub_tag', 'forWeb', '-clb', 'True']
        proc = Po(cmd,stdout=Pi, stderr=Pi)
        procs.append(proc)
for proc in procs:
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
    out_fn_dict[exp] = Ldir['LOo'] / 'tracks2' / Ldir['gtagex'] / (exp + '_' + three_tag + '_forWeb') / ('release_' + dsr + '.nc')
    out_json_dict0[exp] = Ldir['LOo'] / 'tracks2' / Ldir['gtagex'] / (exp + '_' + three_tag + '_forWeb') / ('tracks.json')
    out_json_dict1[exp] = Ldir['LOo'] / 'tracks2' / Ldir['gtagex'] / (exp + '_' + three_tag + '_forWeb') / ('times.json')
result = 'SUCCESS'
for exp in exp_list:
    if not out_fn_dict[exp].is_file():
        result = 'FAIL'
# END OF TRACKER JOBS

# CONVERT TO JSON AND SCP TO HOMER
for exp in exp_list:

    in_fn = out_fn_dict[exp]
    out_fn0 = out_json_dict0[exp]
    out_fn1 = out_json_dict1[exp]

    ds = xr.open_dataset(in_fn)
    # packed time, particle
    lon = ds['lon'].values
    lat = ds['lat'].values
    NT, NP = lon.shape

    # Make a time vector (Note the Time is hours from the start of the release)
    dt0 = datetime.strptime(dsr, Ldir['ds_fmt'])
    Time = ds.Time.values
    dt_list = []
    for h in Time:
        dt_list.append(dt0 + timedelta(days=h/24))

    # Create and save jsons.
    xy = []
    for pp in range(NP):
        xy.append({'x': [('%0.3f' % (item)) for item in lon[:,pp]], 'y': [('%0.3f' % (item)) for item in lat[:,pp]]})
    json.dump(xy, open(out_fn0, 'w'))

    tt_list = []
    for tt in range(NT):
        tt_list.append(datetime.strftime(dt_list[tt],'%Y-%m-%d %H:%M'))
    tt = [{'t': tt_list}]
    json.dump(tt, open(out_fn1, 'w'))
    
    if Ldir['testing']== False:
        
        # send to homer
        cmd2 = ['scp',str(out_json_dict0[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/wgh/tracks.json']
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving tracks to homer ' + out_json_dict[exp].name)
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
            result = 'FAIL'

        cmd2 = ['scp',str(out_json_dict1[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/wgh/times.json']
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving times to homer ' + out_json_dict[exp].name)
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
