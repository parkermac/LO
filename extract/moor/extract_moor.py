"""
This is code for doing mooring extractions.

Test on mac in ipython:

run extract_moor.py -g cas6 -t v3 -x lo8b -ro 1 -0 2019.07.04 -1 2019.07.06

"""

from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun

Ldir = exfun.intro() # this handles the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import Lfun
from time import time
from datetime import timedelta
import subprocess
import netCDF4 as nc
import numpy as np

lon = 10
lat = 10

out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor'
Lfun.make_dir(out_dir)

in_dir0 = Ldir['roms_out'] / Ldir['gtagex']

dt0 = datetime.strptime(Ldir['ds0'], Ldir['ds_fmt'])
dt1 = datetime.strptime(Ldir['ds1'], Ldir['ds_fmt'])

dt = dt0
counter = 0
while dt <= dt1:
    
    tt0 = time()
    
    # proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = proc.communicate()
    
    
    f_string = 'f' + dt.strftime(Ldir['ds_fmt'])
    in_dir = in_dir0 / f_string
    out_fn = out_dir / ('moor_temp_' + str(counter) + '.nc')
    cmd_list = ['ls', str(in_dir)]
    p1 = subprocess.Popen(cmd_list, stdout=subprocess.PIPE)
    cmd_list = ['grep','ocean_his']
    p2 = subprocess.Popen(cmd_list, stdin=p1.stdout, stdout=subprocess.PIPE)
    cmd_list = ['grep','-v','0001']
    p3 = subprocess.Popen(cmd_list, stdin=p2.stdout, stdout=subprocess.PIPE)
    cmd_list = ['ncrcat','-p', str(in_dir),
        '-d', 'xi_rho,'+str(lon), '-d', 'eta_rho,'+str(lat),
        '-d', 'xi_u,'+str(lon), '-d', 'eta_u,'+str(lat),
        '-d', 'xi_v,'+str(lon), '-d', 'eta_v,'+str(lat),
        '-O', str(out_fn)]
    proc = subprocess.Popen(cmd_list, stdin=p3.stdout, stdout=subprocess.PIPE)
    proc.communicate()
    
    cmd_list = ['ls', str(out_dir)]
    pp1 = subprocess.Popen(cmd_list, stdout=subprocess.PIPE)
    cmd_list = ['grep','moor_temp']
    pp2 = subprocess.Popen(cmd_list, stdin=pp1.stdout, stdout=subprocess.PIPE)
    cmd_list = ['ncrcat','-p', str(out_dir), '-O','all.nc']
    proc = subprocess.Popen(cmd_list, stdin=pp2.stdout, stdout=subprocess.PIPE)
    proc.communicate()
    
    counter += 1
    dt += timedelta(days=1)
    print('took %0.2f sec' % (time()-tt0))
    sys.stdout.flush()
    
a = nc.Dataset(out_dir / 'all.nc')
ot = a['ocean_time'][:]
print(np.diff(ot))

# ls | grep ocean_his | ncrcat -d xi_rho,10 -d eta_rho,10 -d xi_u,10 -d eta_u,10 -d xi_v,10 -d eta_v,10 -O moor.nc
# -------------------------------------------------------

# test for success 
if True:
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
