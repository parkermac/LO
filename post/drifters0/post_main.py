"""
This is the main program for doing a couple of tracker runs, converting them
to JSON, and pushing them to my website for some interactive js.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job drifters0 -test True

Run for real on apogee:

run post_main.py -gtx cas6_v0_u0kb -ro 0 -d ['today's datestring] -job drifters0

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
if Ldir['testing']:
    dtt = 1
    dsr = Ldir['date_string']
else:
    dtt = 14
    dtr00 = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    dtr0 = dtr00 - timedelta(days = 11)
    dsr = dtr0.strftime('%Y.%m.%d')

# RUN TRACKER JOBS - parallelize using subprocess
tt0 = time()
procs = []
exp_list = ['full', 'PS']
for exp in exp_list:
        sleep(1) # cludge: needed so that calls to Lstart() don't collide while writing lo_info.csv
        # BUT this is obsolete since we no longer write lo_info.csv 2022.11.08.
        cmd = ['python', str(Ldir['LO']) + '/tracker/tracker.py',
            '-gtx', Ldir['gtagex'], '-ro', str(Ldir['roms_out_num']),
            '-d', dsr, '-dtt', str(dtt),
            '-exp', exp,
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
out_json_dict = {}
out_restructured_json_dict = {}
for exp in exp_list:
    out_fn_dict[exp] = Ldir['LOo'] / 'tracks' / (exp + '_surf_forWeb') / ('release_' + dsr + '.nc')
    out_json_dict[exp] = Ldir['LOo'] / 'tracks' / (exp + '_surf_forWeb') / ('tracks_' + exp + '.json')
    out_restructured_json_dict[exp] = Ldir['LOo'] / 'tracks' / (exp + '_surf_forWeb') / ('tracks_' + exp + '_restructured.json')
result = 'success'
for exp in exp_list:
    if not out_fn_dict[exp].is_file():
        result = 'fail'
# END OF TRACKER JOBS

# CONVERT TO JSON AND SCP TO HOMER
for exp in exp_list:
    if exp == 'full':
        skp = 3
    else:
        skp = 2
    fn = out_fn_dict[exp]
    ds = xr.open_dataset(fn)
    # packed time, particle
    x = ds['lon'].values
    y = ds['lat'].values
    NT, NP = x.shape
    
    # original json
    xy = []
    for pp in range(NP):
        xy.append({'x': list(x[::skp,pp]), 'y': list(y[::skp,pp])})
    json.dump(xy, open(out_json_dict[exp], 'w'))
    
    # restructured json
    rxy = []
    iit_list = list(range(NT))[::skp]
    for iip in range(NP):
        i_newtime = 0
        for iit in iit_list:
            rxy.append({'track':iip, 'point':i_newtime, 'x':round(x[iit,iip],3), 'y':round(y[iit,iip],3)})
            i_newtime += 1
    json.dump(rxy, open(out_restructured_json_dict[exp], 'w'))
    
    if Ldir['testing']== False:
        
        # send original to homer
        cmd2 = ['scp',str(out_json_dict[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks/'+ out_json_dict[exp].name]
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving original json to homer ' + out_json_dict[exp].name)
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
            result = 'fail'
            
        # send restructured to homer
        cmd2 = ['scp',str(out_json_dict[exp]),
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks/'+ out_restructured_json_dict[exp].name]
        proc = Po(cmd2,stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('WARNING: problem moving restructured json to homer ' + out_restructured_json_dict[exp].name)
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
            result = 'fail'
        
    else:
        print('Skipped sending files to homer')
# END CONVERT AND SCP
# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
