"""
This is the main program for making the daymovies and pushing them to homer.

Testing on mac:
run post_main.py -gtx cas6_v3_lo8b -ro 2 -d 2019.07.04 -job daymovie0 -test True

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
from time import time, sleep
import shutil
from lo_tools import Lfun

print((' Creating daymovies for ' + Ldir['date_string'] + ' ').center(60,'='))

ds00 = Ldir['date_string']
dt00 = datetime.strptime(ds00,'%Y.%m.%d')

# if Ldir['testing'] == True:
#     moviename_list = ['P1_full_oxygen_bot']
# else:
moviename_list = ['P1_full_salt_top', 'P1_full_oxygen_bot', 'P1_nshelf_oxygen_bot',
                    'P1_PS_temp_top', 'P1_PS_speed_top',
                    'P1_willapa_ARAG_top', 'P1_willapa_ARAG_bot',
                    'Phab_full_salt_top']

tt0 = time()
result = 'success'
procs = []
for moviename in moviename_list:
    print(moviename.center(60,'-'))
    sys.stdout.flush()
    
    (pt,dom,vn,BOT) = moviename.split('_')
    # set initial assumptions
    tracks = 'False'
    emask = 'False'
    avl = 'True'
    bot = 'False'
    mov = 'True'
    lt = 'hourly'
    # modify assumptions
    if BOT == 'bot':
        bot = 'True'
    if vn in ['oxygen', 'ARAG', 'speed']:
        avl = 'False'
    if (vn in ['oxygen', 'ARAG']) and (dom == 'full'):
        emask = 'True'
    if vn == 'salt':
        tracks = 'True'
    
    if moviename == 'Phab_full_salt_top':
        ttag = 'hab'
    else:
        ttag = 'base'
    
    if moviename == 'P1_nshelf_oxygen_bot':
        # start this one a couple day s earlier
        dt0 = dt00 - timedelta(days=2) # start earlier
        dt1 = dt00 + timedelta(days=Ldir['forecast_days']-1)
    else:
        dt0 = dt00
        dt1 = dt0 + timedelta(days=Ldir['forecast_days']-1)
    ds0 = dt0.strftime('%Y.%m.%d')
    ds1 = dt1.strftime('%Y.%m.%d')
    
    if Ldir['testing']:
        ds1 = ds0
    #     lt = 'snapshot'
    #     mov = 'False'
        
    cmd = ['python', str(Ldir['LO'] / 'daymovie' / 'dm_plot.py'),
        '-gtx', str(Ldir['gtagex']), '-ro', str(Ldir['roms_out_num']),
        '-ds0', ds0, '-ds1', ds1, '-lt', lt, '-mov', mov, '-pt', pt,
        '-dom', dom, '-vn', vn, '-tracks', tracks, '-emask', emask, '-ttag', ttag,
        '-avl', avl, '-bot', bot]
    
    proc = Po(cmd,stdout=Pi, stderr=Pi)
    procs.append(proc)
    sleep(2)

for proc in procs:
    stdout, stderr = proc.communicate()
    if True:
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
    
print('Time to run all jobs = %0.1f sec' % (time() - tt0))
sys.stdout.flush()

out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
# out_dir was already created by post_argfun.py, or by driver_post.py

result = 'success'
for moviename in moviename_list:
    input_filename = Ldir['LOo'] / 'daymovie' / Ldir['gtagex'] / moviename / 'movie.mp4'
    output_filename = moviename + '.mp4'

    if Ldir['testing'] == False:
        # send file to homer (only works from apogee)
        print(' - copying '+output_filename+' to homer')
        sys.stdout.flush()
    
        try:
            cmd_list = ['scp',input_filename,
                'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast/'+output_filename]
            proc = Po(cmd_list,stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                print(' sdtout '.center(60,'-'))
                print(stdout.decode())
            if len(stderr) > 0:
                print(' stderr '.center(60,'-'))
                print(stderr.decode())
            sys.stdout.flush()
        except Exception as e:
            print(' error saving movie to homer')
            print(e)
            result = 'fail'
    else:
        print(' - testing: did not copy files to homer')

    # and save a local copy
    try:
        shutil.copyfile(input_filename, out_dir / output_filename)
    except Exception as e:
        print(' error saving local copy of movie')
        print(e)
        
    sleep(2)
        
# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
