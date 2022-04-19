"""
Code to extract a box-like region, typically for another modeler to use
as a boundary contition.  In cases where it gets velocity in addition to
the rho-grid variables the grid limits mimic the standard ROMS organization,
with the outermost corners being on the rho-grid.

Job definitions are in LO_user/extract/box/job_definitions.py

Testing:
run extract_box -gtx cas6_v3_lo8b -job sequim0 -test True

same but with all flags:
run extract_box -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt daily -job sequim0 -test True

this command replicates what post/surface0 does
run extract_box -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -lt hourly -job surface0 -uv_to_rho True -surf True
or
python extract_box.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -lt hourly -job surface0 -uv_to_rho True -surf True

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
# set this to True to interpolate all u, and v fields to the rho-grid
parser.add_argument('-uv_to_rho', default=False, type=Lfun.boolean_string)
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
    bname = Ldir['job'] + '_surf_' + Ldir['ds0'] + '_' + Ldir['ds1']
    box_fn = out_dir / (bname + '.nc')
elif Ldir['bot']:
    bname = Ldir['job'] + '_bot_' + Ldir['ds0'] + '_' + Ldir['ds1']
    box_fn = out_dir / (bname + '.nc')
else:
    bname = Ldir['job'] + '_' + Ldir['ds0'] + '_' + Ldir['ds1']
    box_fn = out_dir / (bname + '.nc')
box_fn.unlink(missing_ok=True)

# name the temp dir to accumulate individual extractions
temp_dir = out_dir / ('temp_' + bname)
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

# get the job_definitions module, looking first in LO_user
pth = Ldir['LO'] / 'extract' / 'box'
upth = Ldir['LOu'] / 'extract' / 'box'
if (upth / 'job_definitions.py').is_file():
    print('Importing job_definitions from LO_user')
    job_definitions = Lfun.module_from_file('job_definitions', upth / 'job_definitions.py')
else:
    print('Importing job_definitions from LO')
    job_definitions = Lfun.module_from_file('job_definitions', pth / 'job_definitions.py')

aa, vn_list = job_definitions.get_box(Ldir['job'], Lon, Lat)
lon0, lon1, lat0, lat1 = aa
ilon0, ilat0 = check_bounds(lon0, lat0)
ilon1, ilat1 = check_bounds(lon1, lat1)

# NOTE: ncks indexing is zero-based but is INCLUSIVE of the last point.
# NOTE: ncks extractions retain singleton dimensions

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
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            proc.communicate()
        # make sure everyone is finished before continuing
        proc_list = []
    ii += 1

# Ensure that all days have the same fill value.  This was required for cas6_v3_lo8b
# when passing from 2021.10.31 to 2021.11.01 because they had inconsistent fill values,
# which leaks through the ncrcat call below.
tt1 = time()
enc_dict = {'_FillValue':1e20}
vn_List = vn_list.split(',')
Enc_dict = {vn:enc_dict for vn in vn_List}
for out_fn in list(temp_dir.glob('box_*.nc')):
    ds = xr.load_dataset(out_fn) # need to load, not open, for overwrite
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()
print(' - Time for adding fill value = %0.2f sec' % (time()- tt1))

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
print('Time for initial extraction = %0.2f sec' % (time()- tt0))

# add z variables
if (Ldir['surf']==False) and (Ldir['bot']==False):
    tt0 = time()
    ds = xr.load_dataset(box_fn) # have to load in order to add new variables
    NT, N, NR, NC = ds.salt.shape
    ds['z_rho'] = (('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N, NR, NC)))
    ds['z_w'] = (('ocean_time', 's_w', 'eta_rho', 'xi_rho'), np.nan*np.ones((NT, N+1, NR, NC)))
    ds.z_rho.attrs = {'units':'m', 'long_name': 'vertical position on s_rho grid, positive up'}
    ds.z_w.attrs = {'units':'m', 'long_name': 'vertical position on s_w grid, positive up'}
    for ii in range(NT):
        h = ds.h.values
        zeta = ds.zeta[ii,:,:].values
        z_rho, z_w = zrfun.get_z(h, zeta, S)
        ds['z_rho'][ii,:,:,:] = z_rho
        ds['z_w'][ii,:,:,:] = z_w
    ds.to_netcdf(box_fn)
    ds.close()
    print('Time to add z variables = %0.2f sec' % (time()- tt0))

if Ldir['uv_to_rho']:
    # interpolate anything on the u and v grids to the rho grid, assuming
    # zero values where masked, and leaving a masked ring around the outermost edge
    tt0 = time()
    ds = xr.load_dataset(box_fn) # have to load in order to add new variables
    Maskr = ds.mask_rho.values == 1 # True over water
    NR, NC = Maskr.shape
    for vn in ds.data_vars:
        if ('xi_u' in ds[vn].dims) and ('ocean_time' in ds[vn].dims):
            if len(ds[vn].dims) == 4:
                uu = ds[vn].values
                NT, N, NRu, NCu = uu.shape
                uu[np.isnan(uu)] = 0
                UU = (uu[:,:,1:-1,1:]+uu[:,:,1:-1,:-1])/2
                uuu = np.nan * np.ones((NT, N, NR, NC))
                uuu[:,:,1:-1,1:-1] = UU
                Maskr3 = np.tile(Maskr.reshape(1,1,NR,NC),[NT,N,1,1])
                uuu[~Maskr3] = np.nan
                ds.update({vn:(('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), uuu)})
            elif len(ds[vn].dims) == 3:
                uu = ds[vn].values
                NT, NRu, NCu = uu.shape
                uu[np.isnan(uu)] = 0
                UU = (uu[:,1:-1,1:]+uu[:,1:-1,:-1])/2
                uuu = np.nan * np.ones((NT, NR, NC))
                uuu[:,1:-1,1:-1] = UU
                Maskr3 = np.tile(Maskr.reshape(1,NR,NC),[NT,1,1])
                uuu[~Maskr3] = np.nan
                ds.update({vn:(('ocean_time', 'eta_rho', 'xi_rho'), uuu)})
        elif ('xi_v' in ds[vn].dims) and ('ocean_time' in ds[vn].dims):
            if len(ds[vn].dims) == 4:
                vv = ds[vn].values
                NT, N, NRv, NCv = vv.shape
                vv[np.isnan(vv)] = 0
                VV = (vv[:,:,1:,1:-1]+vv[:,:,:-1,1:-1])/2
                vvv = np.nan * np.ones((NT, N, NR, NC))
                vvv[:,:,1:-1,1:-1] = VV
                Maskr3 = np.tile(Maskr.reshape(1,1,NR,NC),[NT,N,1,1])
                vvv[~Maskr3] = np.nan
                ds.update({vn:(('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), vvv)})
            elif len(ds[vn].dims) == 3:
                vv = ds[vn].values
                NT, NRv, NCv = vv.shape
                vv[np.isnan(vv)] = 0
                VV = (vv[:,1:,1:-1]+vv[:,:-1,1:-1])/2
                vvv = np.nan * np.ones((NT, NR, NC))
                vvv[:,1:-1,1:-1] = VV
                Maskr3 = np.tile(Maskr.reshape(1,NR,NC),[NT,1,1])
                vvv[~Maskr3] = np.nan
                ds.update({vn:(('ocean_time', 'eta_rho', 'xi_rho'), vvv)})
    ds.to_netcdf(box_fn)
    ds.close()
    print('Time to interpolate uv variables to rho grid = %0.2f sec' % (time()- tt0))

# squeeze and compress the resulting file
tt0 = time()
ds = xr.load_dataset(box_fn)
ds = ds.squeeze() # remove singleton dimensions
enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
Enc_dict = {vn:enc_dict for vn in ds.data_vars if 'ocean_time' in ds[vn].dims}
ds.to_netcdf(box_fn, encoding=Enc_dict)
ds.close()
print('Time to compress = %0.2f sec' % (time()- tt0))

# clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()

print('Size of full rho-grid = %s' % (str(G['lon_rho'].shape)))
print(' Contents of extracted box file: '.center(60,'-'))
# check on the results
ds = xr.open_dataset(box_fn)
for vn in ds.data_vars:
    print('%s %s max/min = %0.4f/%0.4f' % (vn, str(ds[vn].shape), ds[vn].max(), ds[vn].min()))
ds.close()

print('\nPath to file:\n%s' % (str(box_fn)))


