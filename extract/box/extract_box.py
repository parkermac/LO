"""
Code to extract a box-like region, typically for another modeler to use
as a boundary contition.

Testing:
run extract_box -gtx cas6_v3_lo8b -job yang_sequim -test True

same but with all flags:
run extract_box -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -lt hourly -job yang_sequim -test True

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
    Ldir['ds1'] = '2019.07.04'
    Ldir['list_type'] = 'hourly'
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

def get_box(job):
    vn_list = 'salt,temp,zeta,h,u,v' # default list
    # specific jobs
    if job == 'yang_sequim':
        aa = [-123.15120787, -122.89090010, 48.07302111, 48.19978336]
    elif job == 'PS':
        # 3 MB per save (26 GB/year for hourly)
        aa = [-123.5, -122.05, 47, 49]
        Ldir['surf'] = True # override to be sure
    elif job == 'garrison':
        aa = [-129.9, -122.05, 42.1, 51.9]
        Ldir['bot'] = True # override to be sure
        vn_list = 'salt,temp,oxygen,h'
    return aa[0], aa[1], aa[2], aa[3], vn_list
    
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
lon0, lon1, lat0, lat1, vn_list = get_box(Ldir['job'])
ilon0, ilat0 = check_bounds(lon0, lat0)
ilon1, ilat1 = check_bounds(lon1, lat1)

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
        '-d', 'xi_u,'+str(ilon0-1)+','+str(ilon1), '-d', 'eta_u,'+str(ilat0)+','+str(ilat1),
        '-d', 'xi_v,'+str(ilon0)+','+str(ilon1), '-d', 'eta_v,'+str(ilat0-1)+','+str(ilat1)]
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
    ds['z_rho'] = 0*ds.temp
    ds.z_rho.attrs['units'] = 'm'
    ds.z_rho.attrs['long name'] = 'vertical position on s_rho grid, positive up'
    NT = len(ds.ocean_time.values)
    for ii in range(NT):
        h = ds.h.values
        zeta = ds.zeta[ii,:,:].values
        z_rho = zrfun.get_z(h, zeta, S, only_rho=True)
        ds['z_rho'][ii,:,:,:] = z_rho
    ds.to_netcdf(box_fn)
    ds.close()
    print('time to add z variables = %0.2f sec' % (time()- tt0))

# clean up
if True: #Ldir['testing'] == False:
    Lfun.make_dir(temp_dir, clean=True)
    temp_dir.rmdir()

print('\nContents of extracted box file:')
# check on the results
ds = xr.open_dataset(box_fn)
for vn in ds.data_vars:
    print('%s (%s) max/min = %0.4f/%0.4f' % (vn, str(ds[vn].shape), ds[vn].max(), ds[vn].min()))
ds.close()

# ** this should be moved to a separate plotting program **
if Ldir['testing'] and Ldir['lo_env'] == 'pm_mac':
    # make a plot to look at things
    # (currently optimized for yang_sequim)
    import matplotlib.pyplot as plt
    import plotting_functions as pfun
    ds = xr.open_dataset(box_fn)
    s0 = ds.salt[0,-1,:,:].values
    u0 = ds.u[0,-1,:,:].values
    v0 = ds.v[0,-1,:,:].values
    rmask = ~np.isnan(s0) # True on water
    umask = ~np.isnan(u0)
    vmask = ~np.isnan(v0)
    
    plt.close('all')
    pfun.start_plot(figsize=(22,9))
    fig = plt.figure()
    
    ax = fig.add_subplot(121)
    alpha=.2
    ms = 3
    ax.plot(ds.lon_rho.values, ds.lat_rho.values,'ok',alpha=alpha, ms=ms)
    ax.plot(ds.lon_u.values, ds.lat_u.values,'>g',alpha=alpha, ms=ms)
    ax.plot(ds.lon_v.values, ds.lat_v.values,'^r',alpha=alpha, ms=ms)
    ax.plot(ds.lon_rho.values[rmask], ds.lat_rho.values[rmask],'ok', ms=ms)
    ax.plot(ds.lon_u.values[umask], ds.lat_u.values[umask],'>g', ms=ms)
    ax.plot(ds.lon_v.values[vmask], ds.lat_v.values[vmask],'^r', ms=ms)
    pfun.add_coast(ax, color='b', linewidth=2)
    pfun.dar(ax)
    pad = .02
    ax.axis([lon0-pad, lon1+pad, lat0-pad, lat1+pad])
    
    # make psi_grid coordinates
    lon_psi, lat_psi = np.meshgrid(ds.lon_u.values[0,:], ds.lat_v.values[:,0])
    # make interpolated u and v
    uu = u0.copy(); vv = v0.copy()
    uu[np.isnan(uu)] = 0
    vv[np.isnan(vv)] = 0
    UU = (uu[:,1:]+uu[:,:-1])/2
    VV = (vv[:1,:] + vv[:-1,:])/2
    UU[np.isnan(s0)] = np.nan
    VV[np.isnan(s0)] = np.nan
    
    ax = fig.add_subplot(122)
    alpha=.2
    ax.pcolormesh(lon_psi, lat_psi, s0, vmin=30, vmax=32, cmap='Spectral_r')
    ax.quiver(ds.lon_rho.values, ds.lat_rho.values, UU, VV)
    pfun.add_coast(ax, color='b', linewidth=2)
    pfun.dar(ax)
    pad = .02
    ax.axis([lon0-pad, lon1+pad, lat0-pad, lat1+pad])
    
    plt.show()
    pfun.end_plot()
    
    ds.close()
