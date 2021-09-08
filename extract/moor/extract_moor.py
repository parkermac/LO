"""
This is code for doing mooring extractions.

Test on mac in ipython:
run extract_moor -gtx cas6_v3_lo8b -test True

The same job would be run with flags as:
python extract_moor.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt hourly -sn TEST_0 -lon ' -125' -lat 47 -get_all True > em.log &
NOTE: the quotes and space are required to feed it a negative longitude.

The performance on this is excellent, taking about 24 minutes for a year of hourly records
on perigee with cas6_v3_lo8b and all flags True, IF we use Nproc = 100.

Using the more standard Nproc = 10, it takes 1.2 hours for a year on perigee.

"""

# imports
import sys
from lo_tools import Lfun, zrfun, zfun
import argparse
from time import time
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import numpy as np
import xarray as xr

# commandd line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 1 = Ldir['roms_out1'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly or daily
# select mooring location and name
parser.add_argument('-lon', type=float) # longitude
parser.add_argument('-lat', type=float) # latitude
parser.add_argument('-sn', type=str)    # station name
# select categories of variables to extract (defined in extract_moor.py)
parser.add_argument('-get_tsa', type=zfun.boolean_string, default=False)
parser.add_argument('-get_vel', type=zfun.boolean_string, default=False)
parser.add_argument('-get_bio', type=zfun.boolean_string, default=False)
parser.add_argument('-get_surfbot', type=zfun.boolean_string, default=False)
# OR select all of them
parser.add_argument('-get_all', type=zfun.boolean_string, default=False)
# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)    
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided (other)
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument to forcing_argfun.intro(): ' + a)
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
    Ldir['list_type'] = 'hourly'
    Ldir['lon'] = -125
    Ldir['lat'] = 47
    Ldir['sn'] = 'TEST_0'
    Ldir['get_all'] = True
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
# set variable list flags
if Ldir['get_all']:
    Ldir['get_tsa'] = True
    Ldir['get_vel'] = True
    Ldir['get_bio'] = True
    Ldir['get_surfbot'] = True

# do the extraction
tt00 = time()

# set output location
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor'
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor' / ('temp_' + Ldir['sn'])
Lfun.make_dir(out_dir)
Lfun.make_dir(temp_dir, clean=True)
moor_fn = out_dir / (Ldir['sn'] + '_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
moor_fn.unlink(missing_ok=True)
print(moor_fn)

# get indices for extraction
in_dir0 = Ldir['roms_out'] / Ldir['gtagex']
lon = Ldir['lon']
lat = Ldir['lat']
G, S, T = zrfun.get_basic_info(in_dir0 / ('f' + Ldir['ds0']) / 'ocean_his_0001.nc')
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]    
# error checking
if (lon < Lon[0]) or (lon > Lon[-1]):
    print('ERROR: lon out of bounds ' + moor_fn.name)
    sys.exit()
if (lat < Lat[0]) or (lat > Lat[-1]):
    print('ERROR: lat out of bounds ' + moor_fn.name)
    sys.exit()
# get indices
ilon = zfun.find_nearest_ind(Lon, lon)
ilat = zfun.find_nearest_ind(Lat, lat)
# more error checking
if G['mask_rho'][ilat,ilon] == False:
    print('ERROR: rho point on land mask ' + moor_fn.name)
    sys.exit()
if Ldir['get_vel'] or Ldir['get_surfbot']:
    if G['mask_u'][ilat,ilon] == False:
        print('ERROR: u point on land mask ' + moor_fn.name)
        sys.exit()
    if G['mask_v'][ilat,ilon] == False:
        print('ERROR: v point on land mask ' + moor_fn.name)
        sys.exit()
        
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

vn_list = 'h,zeta'
if Ldir['get_tsa']:
    vn_list += ',salt,temp,AKs,AKv'
if Ldir['get_vel']:
    vn_list += ',u,v,w,ubar,vbar'
if Ldir['get_bio']:
    vn_list += ',NO3,phytoplankton,zooplankton,detritus,Ldetritus,oxygen,alkalinity,TIC'
if Ldir['get_surfbot']:
    vn_list += ',Pair,Uwind,Vwind,shflux,ssflux,latent,sensible,lwrad,swrad,sustr,svstr,bustr,bvstr'
    
proc_list = []
N = len(fn_list)
print('Times to extract =  %d' % (N))
for ii in range(N):
    fn = fn_list[ii]
    # extract one day at a time using ncks
    count_str = ('000000' + str(ii))[-6:]
    out_fn = temp_dir / ('moor_temp_' + count_str + '.nc')
    cmd_list1 = ['ncks',
        '-v', vn_list,
        '-d', 'xi_rho,'+str(ilon), '-d', 'eta_rho,'+str(ilat)]
    if Ldir['get_vel'] or Ldir['get_surfbot']:
        cmd_list1 += ['-d', 'xi_u,'+str(ilon), '-d', 'eta_u,'+str(ilat),
            '-d', 'xi_v,'+str(ilon), '-d', 'eta_v,'+str(ilat)]
    cmd_list1 += ['-O', str(fn), str(out_fn)]
    proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    
    if (np.mod(ii,100) == 0): # 100
        print(str(ii), end=', ')
        sys.stdout.flush()
        if (np.mod(ii,1000) == 0) and (ii > 0): # 1000
            print(str(ii))
            sys.stdout.flush()
    elif (ii == N-1):
        print(str(ii))
        sys.stdout.flush()
    
    # Nproc controls how many ncks subprocesses we allow to stack up
    # before we require them all to finish.  It appears to work even
    # with Nproc = 100 on perigee, although this may slow other jobs.
    # boiler seemed to prefer a lower number, like 10.
    Nproc = Ldir['Nproc']
    if (np.mod(ii,Nproc) == 0) or (ii == N-1):
        for proc in proc_list:
            proc.communicate()
        # make sure everyone is finished before continuing
        proc_list = []

# concatenate the day records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','moor_temp'], stdin=pp1.stdout, stdout=Pi)
cmd_list = ['ncrcat','-p', str(temp_dir), '-O',str(moor_fn)]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi)
proc.communicate()

# add z_coordinates to the file using xarray
ds = xr.load_dataset(moor_fn)
ds = ds.squeeze() # remove singleton dimensions
zeta = ds.zeta.values
NT = len(zeta)
hh = ds.h.values * np.ones(NT)
z_rho, z_w = zrfun.get_z(hh, zeta, S)
# the returned z arrays have vertical position first, so we 
# transporse to put time first for the mooring, to be consistent with
# all other variables
ds['z_rho'] = (('ocean_time', 's_rho'), np.transpose(z_rho.data))
ds['z_w'] = (('ocean_time', 's_w'), np.transpose(z_w.data))
ds.z_rho.attrs['units'] = 'm'
ds.z_w.attrs['units'] = 'm'
ds.z_rho.attrs['long name'] = 'vertical position on s_rho grid, positive up'
ds.z_w.attrs['long name'] = 'vertical position on s_w grid, positive up'
# add units to salt
ds.salt.attrs['units'] = 'g kg-1'
# update the time long name
ds.ocean_time.attrs['long_name'] = 'Time [UTC]'
# update format attribute
ds.attrs['format'] = 'netCDF-4'
# and save to NetCDF (default is netCDF-4, and to overwrite any existing file)
ds.to_netcdf(moor_fn)
ds.close()
    
# clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()

print('Total Elapsed time was %0.2f sec' % (time()-tt00))
