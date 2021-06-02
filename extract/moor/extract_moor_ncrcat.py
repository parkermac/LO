"""
This is code for doing mooring extractions.

Test on mac in ipython:

run extract_moor.py -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06

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
import zrfun, zfun
from time import time
from datetime import timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi

# set output location
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor'
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor' / ('temp_' + Ldir['sn'])
Lfun.make_dir(out_dir)
Lfun.make_dir(temp_dir, clean=True)
moor_fn = out_dir / (Ldir['sn'] + '.nc')
moor_fn.unlink(missing_ok=True)

# get indices for extraction
in_dir0 = Ldir['roms_out'] / Ldir['gtagex']
lon = Ldir['lon']
lat = Ldir['lat']
G = zrfun.get_basic_info(in_dir0 / ('f' + Ldir['ds0']) / 'ocean_his_0001.nc', only_G=True)
ilon = zfun.find_nearest_ind(G['lon_rho'][0,:], lon)
ilat = zfun.find_nearest_ind(G['lat_rho'][:,0], lat)
# NOTE: we should also check that this is not masked on any grid

# do the extraction
dt0 = datetime.strptime(Ldir['ds0'], Ldir['ds_fmt'])
dt1 = datetime.strptime(Ldir['ds1'], Ldir['ds_fmt'])
dt = dt0
counter = 0
while dt <= dt1:
    tt0 = time()
    
    # extract one day at a time using ncrcat
    f_string = 'f' + dt.strftime(Ldir['ds_fmt'])
    in_dir = in_dir0 / f_string
    count_str = ('0000' + str(counter))[-4:]
    out_fn = temp_dir / ('moor_temp_' + count_str + '.nc')
    p1 = Po(['ls', str(in_dir)], stdout=Pi)
    p2 = Po(['grep','ocean_his'], stdin=p1.stdout, stdout=Pi)
    cmd_list = ['ncrcat','-p', str(in_dir),
        '-d', 'xi_rho,'+str(ilon), '-d', 'eta_rho,'+str(ilat),
        '-d', 'xi_u,'+str(ilon), '-d', 'eta_u,'+str(ilat),
        '-d', 'xi_v,'+str(ilon), '-d', 'eta_v,'+str(ilat),
        '-O', str(out_fn)]
    if dt == dt0:
        # on the first day do all hours
        proc = Po(cmd_list, stdin=p2.stdout, stdout=Pi)
    else:
        # on subsequent days skip the (repeated) hour 0
        p3 = Po(['grep','-v','0001'], stdin=p2.stdout, stdout=Pi)
        proc = Po(cmd_list, stdin=p3.stdout, stdout=Pi)
    proc.communicate()
    
    counter += 1
    dt += timedelta(days=1)
    
    print('took %0.2f sec' % (time()-tt0))
    sys.stdout.flush()
    
# concatenate the day records into one file
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','moor_temp'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O',str(moor_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi)
proc.communicate()
    
# clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()

# test for success 
if moor_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
