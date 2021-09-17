"""
Code to extract a box-like region, typically for another modeler to use
as a boundary contition.  In cases where it gets velocity in addition to
the rho-grid variables the grid limits mimic the standard ROMS organization,
with the outermost corners being on the rho-grid.

Job definitions are in LO_user/extract/box/job_definitions.py

Testing:
run extract_box -gtx cas6_v3_lo8b -job yang_sequim -test True

same but with all flags:
run extract_box -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt daily -job yang_sequim -test True

Performance: this is very fast, takes just a few seconds for three days on boiler (for yang_sequim).
"""

# imports
import sys
import argparse
from lo_tools import Lfun, zfun, zrfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
from time import time
import numpy as np
import xarray as xr

pid = os.getpid()
print(' extract_box '.center(60,'='))
print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly
# select job name
parser.add_argument('-job', type=str) # job name
# these flags get only surface or bottom fields if True
# - cannot have both True -
parser.add_argument('-surf', default=False, type=Lfun.boolean_string)
parser.add_argument('-bot', default=False, type=Lfun.boolean_string)
# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# testing
if Ldir['testing']:
    Ldir['roms_out_num'] = 2
    Ldir['ds0'] = '2019.07.04'
    Ldir['ds1'] = '2019.07.06'
    Ldir['list_type'] = 'daily'
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
# check for input conflicts:
if Ldir['surf'] and Ldir['bot']:
    print('Error: cannot have surf and bot both True.')
    sys.exit()
    
# output location
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'box'
Lfun.make_dir(out_dir)
if Ldir['surf']:
    box_fn = out_dir / (Ldir['job'] + '_surf_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
elif Ldir['bot']:
    box_fn = out_dir / (Ldir['job'] + '_bot_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
else:
    box_fn = out_dir / (Ldir['job'] + '_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
box_fn.unlink(missing_ok=True)

# name the temp dir to accumulate individual extractions
temp_dir = out_dir / ('temp_' + Ldir['job'])
Lfun.make_dir(temp_dir, clean=True)

# get list of files to work on
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])
if Ldir['testing']:
    fn_list = fn_list[:5]
G, S, T = zrfun.get_basic_info(fn_list[0])
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]
    
def check_bounds(lon, lat):
    # error checking
    if (lon < Lon[0]) or (lon > Lon[-1]):
        print('ERROR: lon out of bounds ')
        sys.exit()
    if (lat < Lat[0]) or (lat > Lat[-1]):
        print('ERROR: lat out of bounds ')
        sys.exit()
    # get indices
    ilon = zfun.find_nearest_ind(Lon, lon)
    ilat = zfun.find_nearest_ind(Lat, lat)
    return ilon, ilat

# get the indices and check that they are in the grid
pth = Ldir['LOu'] / 'extract' / 'box'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import job_definitions
from importlib import reload
reload(job_definitions)
aa, vn_list = job_definitions.get_box(Ldir['job'])
lon0, lon1, lat0, lat1 = aa
ilon0, ilat0 = check_bounds(lon0, lat0)
ilon1, ilat1 = check_bounds(lon1, lat1)
# NOTE: ncks indexing is zero-based but is INCLUSIVE of the last point.

# do the extractions
N = len(fn_list)
proc_list = []
tt0 = time()
print('Working on ' + box_fn.name + ' (' + str(N) + ' times)')
for ii in range(N):
    fn = fn_list[ii]
    sys.stdout.flush()
    # extract one day at a time using ncks
    count_str = ('000000' + str(ii))[-6:]
    out_fn = temp_dir / ('box_' + count_str + '.nc')
    cmd_list1 = ['ncks',
        '-v', vn_list,
        '-d', 'xi_rho,'+str(ilon0)+','+str(ilon1), '-d', 'eta_rho,'+str(ilat0)+','+str(ilat1),
        '-d', 'xi_u,'+str(ilon0)+','+str(ilon1-1), '-d', 'eta_u,'+str(ilat0)+','+str(ilat1),
        '-d', 'xi_v,'+str(ilon0)+','+str(ilon1), '-d', 'eta_v,'+str(ilat0)+','+str(ilat1-1)]
    if Ldir['surf']:
        cmd_list1 += ['-d','s_rho,'+str(S['N']-1)]
    elif Ldir['bot']:
        cmd_list1 += ['-d','s_rho,0']
    cmd_list1 += ['-O', str(fn), str(out_fn)]
    proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
    proc_list.append(proc)

    # screen output about progress
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == N-1):
        print(str(ii))
        sys.stdout.flush()
        
    # Nproc controls how many ncks subprocesses we allow to stack up
    # before we require them all to finish.
    if (np.mod(ii,Ldir['Nproc']) == 0) or (ii == N-1):
        for proc in proc_list:
            proc.communicate()
        # make sure everyone is finished before continuing
        proc_list = []
    ii += 1
    
# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','box'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O', str(box_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if Ldir['testing']:
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())
print('Total time = %0.2f sec' % (time()- tt0))

# add z variables
if (Ldir['surf']==False) and (Ldir['bot']==False):
    tt0 = time()
    ds = xr.load_dataset(box_fn) # have to load in order to add new variables
    NT, N, NR, NC = ds.salt.shape
    ds['z_rho'] = (('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N, NR, NC)))
    ds['z_w'] = (('ocean_time', 's_w', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N+1, NR, NC)))
    ds.z_rho.attrs = {'units':'m', 'long_name': 'vertical position on s_rho grid, positive up'}
    ds.z_rho.attrs = {'units':'m', 'long_name': 'vertical position on s_w grid, positive up'}
    for ii in range(NT):
        h = ds.h.values
        zeta = ds.zeta[ii,:,:].values
        z_rho, z_w = zrfun.get_z(h, zeta, S)
        ds['z_rho'][ii,:,:,:] = z_rho
        ds['z_w'][ii,:,:,:] = z_w
    ds.to_netcdf(box_fn)
    ds.close()
    print('time to add z variables = %0.2f sec' % (time()- tt0))

# clean up
if True:
    Lfun.make_dir(temp_dir, clean=True)
    temp_dir.rmdir()

print('Size of full rho-grid = %s' % (str(G['lon_rho'].shape)))
print('\nContents of extracted box file:')
# check on the results
ds = xr.open_dataset(box_fn)
for vn in ds.data_vars:
    print('%s %s max/min = %0.4f/%0.4f' % (vn, str(ds[vn].shape), ds[vn].max(), ds[vn].min()))
ds.close()


