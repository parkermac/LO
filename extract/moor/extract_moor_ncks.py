"""
This is code for doing mooring extractions.

Test on mac in ipython:

run extract_moor_ncks.py -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -get_tsa True -get_vel True -get_bio True

The performance on this is amazing.
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
import numpy as np
import netCDF4 as nc
from multiprocessing.pool import ThreadPool
import multiprocessing
ncpu = multiprocessing.cpu_count()
print('ncpu = %d' % (ncpu))

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
G, S, T = zrfun.get_basic_info(in_dir0 / ('f' + Ldir['ds0']) / 'ocean_his_0001.nc')
ilon = zfun.find_nearest_ind(G['lon_rho'][0,:], lon)
ilat = zfun.find_nearest_ind(G['lat_rho'][:,0], lat)
# NOTE: we should also check that this is not masked on any grid

fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

vn_list = 'h,zeta'
if Ldir['get_tsa']:
    vn_list += ',salt,temp,AKs,AKv'
if Ldir['get_vel']:
    vn_list += ',u,v,w'
if Ldir['get_bio']:
    vn_list += ',NO3,phytoplankton,zooplankton,detritus,Ldetritus,oxygen,alkalinity,TIC'
if Ldir['get_surfbot']:
    pass
    # need to populate surface fields to get
    

proc_list = []    
def my_fun(ii):
    fn = fn_list[ii]
    # extract one day at a time using ncrcat
    count_str = ('0000' + str(ii))[-4:]
    out_fn = temp_dir / ('moor_temp_' + count_str + '.nc')
    cmd_list1 = ['ncks',
        '-v', vn_list,
        '-d', 'xi_rho,'+str(ilon), '-d', 'eta_rho,'+str(ilat)]
    if Ldir['get_vel']:
        cmd_list1 += ['-d', 'xi_u,'+str(ilon), '-d', 'eta_u,'+str(ilat),
            '-d', 'xi_v,'+str(ilon), '-d', 'eta_v,'+str(ilat)]
    cmd_list1 += ['-O', str(fn), str(out_fn)]
    proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    #counter += 1

tp = ThreadPool(None) # defaults to number of processors
tt0 = time()
for ii in range(len(fn_list)):
    tp.apply_async(my_fun, (ii,))
tp.close()
tp.join()
# make sure everyone is finished before continuing
for proc in proc_list:
    proc.communicate()
print('Total days %3d: took %0.2f sec' % (ii/24, time()-tt0))
sys.stdout.flush()
    
# concatenate the day records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','moor_temp'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O',str(moor_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi)
proc.communicate()

# Add z-coordinates to the file
foo = nc.Dataset(moor_fn, 'a')
z_rho, z_w = zrfun.get_z(foo['h'][:], np.array([0.]), S)
vv = foo.createVariable('z_rho', float, ('s_rho',))
vv.long_name = 'vertical position on s_rho grid, positive up, zero at surface'
vv.units = 'm'
vv[:] = z_rho
vv = foo.createVariable('z_w', float, ('s_w',))
vv.long_name = 'vertical position on s_w grid, positive up, zero at surface'
vv.units = 'm'
vv[:] = z_w
# Note that z_rho and z_w do not have singleton dimensions

# Also add units to salt
foo['salt'].units = 'g kg-1'
foo.close()
    
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
