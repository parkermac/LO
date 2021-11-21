"""
This makes the ocn forcing files for a nested run, interpolating to get
the fields from history files for another run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g hc0 -t v0 -r backfill -s continuation -d 2019.07.04 -f ocnN -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
from scipy.spatial import cKDTree
import numpy as np
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import pickle

from lo_tools import Lfun, zfun, zrfun

import Ofun_nc_xarray
# defaults
verbose = True
if Ldir['testing']:
    # verbose = True
    from importlib import reload
    reload(Ofun_nc_xarray)

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# where to find the files to interpolate from
# NOTE: this should be made more general - perhaps handled by the command line arguments.
if 'apogee' in Ldir['lo_env']:
    in_dir = Ldir['roms_out2'] / 'cas6_v3_lo8b' / ('f' + Ldir['date_string'])
    
else:
    in_dir = Ldir['parent'] / 'LiveOcean_roms' / 'output' / 'cas6_v3_lo8b' / ('f' + Ldir['date_string'])

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# list of history files from the original grid to work from (all Path objects)
h_list = sorted(in_dir.glob('ocean_his_*'))
if Ldir['testing']:
    h_list = h_list[:2]

# meant to be replaced by get_one_time.py
if False:
    
    def get_bounds(x_big, y_big, x_small, y_small, pad=3):
        """
        This function takes two pairs of plaid, 2-D, lon, lat arrays:
        - one pair bigger (that we hope to nest inside) and
        - one pair smaller (the grid of the nest)
        and returns the indices to use for making a trimmed version of
        the bigger grid that the smaller grid still fits inside.
        We add "pad" around the edges to make sure things fit comfortably.
        """
        # First: error checking
        if (x_small[0,0] < x_big[0,pad]) or (x_small[0,-1] > x_big[0,-pad]):
            print('ERROR: lon out of bounds ')
            sys.exit()
        if (y_small[0,0] < y_big[pad,0]) or (y_small[-1,0] > y_big[-pad,0]):
            print('ERROR: lat out of bounds ')
            sys.exit()
        # Second: get indices
        ix0 = zfun.find_nearest_ind(x_big[0,:], x_small[0,0]) - pad
        ix1 = zfun.find_nearest_ind(x_big[0,:], x_small[0,-1]) + pad
        iy0 = zfun.find_nearest_ind(y_big[:,0], y_small[0,0]) - pad
        iy1 = zfun.find_nearest_ind(y_big[:,0], y_small[-1,0]) + pad
        return ix0, ix1, iy0, iy1

    tag_list = ['rho', 'u', 'v']

    # the new grid
    ds = xr.open_dataset(Ldir['grid'] / 'grid.nc')

    xx = {}; yy = {}; mm = {}; xynew = {}
    for tag in tag_list:
        xx[tag] = ds['lon_' + tag].values
        yy[tag] = ds['lat_' + tag].values
        mm[tag] = ds['mask_' + tag].values
    
    ds.close()

    if Ldir['start_type'] == 'continuation':
        pad = 20
    elif Ldir['start_type'] == 'new':
        pad = 0
    else:
        print('Error: Unrecognized start_type')
        sys.exit()
    if pad > 0:
        # mask out the inside of the nest fields, since we only use
        # the edges (unless Ldir['start_type']=='new')
        for tag in tag_list:
            mm[tag][pad:-pad, pad:-pad] = 0 # this speeds things up
    
    for tag in tag_list:
        xynew[tag] = np.array((xx[tag][mm[tag]==1],yy[tag][mm[tag]==1])).T

    tt0 = time()
    # Create 2-D search trees for the old grid
    ds = xr.open_dataset(h_list[0])
    N = len(ds.s_rho.values)
    xtrim = {}; ytrim = {}; mtrim = {}; xyT = {}
    ix0 = {}; ix1 = {}; iy0 = {}; iy1 = {}
    for tag in tag_list:
        x = ds['lon_' + tag].values
        y = ds['lat_' + tag].values
        m = ds['mask_' + tag].values # 1=water
        # trim the old grid before making the search tree
        ix0[tag], ix1[tag], iy0[tag], iy1[tag] = get_bounds(x, y, xx[tag], yy[tag])
        xtrim[tag] = x[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
        ytrim[tag] = y[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
        mtrim[tag] = m[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
        xyorig = np.array((xtrim[tag][mtrim[tag]==1],ytrim[tag][mtrim[tag]==1])).T
        xyT[tag] = cKDTree(xyorig)
    ds.close()
    print('Time to make Trees = %0.2f sec' % (time()-tt0))
    sys.stdout.flush()

    # associate variables to process with grids
    vn_dict = {'salt':('rho',3), 'temp':('rho',3), 'zeta':('rho',2),
            'u':('u',3), 'v':('v',3), 'ubar':('u',2), 'vbar':('v',2)}

    # create blank arrays for results
    data_dict = dict()
    for vn in vn_dict.keys():
        tag = vn_dict[vn][0]
        dm = vn_dict[vn][1]
        NR, NC = xx[tag].shape
        if dm == 2:
            data_dict[vn] = np.nan* np.ones((NT,NR,NC))
        elif dm == 3:
            data_dict[vn] = np.nan* np.ones((NT,N,NR,NC))

    # Interpolate to fill all data arrays for new grid.
    tt = 0
    modtime_list = []
    for fn in h_list:
        tt0 = time()
        ds = xr.open_dataset(fn, decode_times=False)
        modtime_list.append(ds.ocean_time.values[0])
        for vn in vn_dict.keys():
            tag = vn_dict[vn][0]
            dm = vn_dict[vn][1]
            if dm == 2:
                vtrim = ds[vn][0,iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]].values
                vv = np.nan * np.ones(xx[tag].shape) 
                vv[mm[tag]==1] = vtrim[mtrim[tag]==1][xyT[tag].query(xynew[tag], workers=-1)[1]]
                # note that "workers" has replaced "n_jobs"
                data_dict[vn][tt, :, :] = vv
            elif dm == 3:
                for nn in range(N):
                    vtrim = ds[vn][0,nn,iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]].values
                    vv = np.nan * np.ones(xx[tag].shape) 
                    vv[mm[tag]==1] = vtrim[mtrim[tag]==1][xyT[tag].query(xynew[tag], workers=-1)[1]]
                    data_dict[vn][tt, nn, :, :] = vv
        print('tt = %d (%0.2f sec)' % (tt, time()-tt0))
        sys.stdout.flush()
        tt += 1
        ds.close()
    data_dict['ocean_time'] = np.array(modtime_list)


# +++++++++++ parallel subprocess +++++++++++++++++++++++++++++++++++
temp_dir = out_dir / 'Data'
Lfun.make_dir(temp_dir, clean=True)
Nproc = 10
NT = len(h_list)
proc_list = []
tt0 = time()
temp_out_fn_list = []
print('Working on ' + Ldir['frc'] + ' (' + str(NT) + ' times)')
for ii in range(NT):
    grid_fn = str(Ldir['grid'] / 'grid.nc')
    his_fn = str(h_list[ii])
    count_str = ('000000' + str(ii))[-6:]
    temp_out_fn = temp_dir / ('data_dict_' + count_str + '.p')
    temp_out_fn_list.append(temp_out_fn)
    this_dir = str(Ldir['LO'] / 'forcing' / Ldir['frc']) + '/'
    cmd_list = ['python', this_dir + 'get_one_time.py', '-grid_fn', grid_fn,
         '-his_fn', his_fn, '-start_type', Ldir['start_type'], '-out_fn', str(temp_out_fn)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
        
    if ((np.mod(ii,Nproc) == 0) and (ii > 0)) or (ii == NT-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            # make sure everyone is finished before continuing
            if True:
                if len(stdout) > 0:
                    print('\n'+stdout.decode())
                if len(stderr) > 0:
                    print('\n'+stderr.decode())
        proc_list = []
    ii += 1
print('Time to run all extractions = %0.1f sec' % (time()-tt0))
sys.stdout.flush()

# concatenate the dicts into one file
data_dict = dict()
ii = 0
for temp_out_fn in temp_out_fn_list:
    dd = pickle.load(open(temp_out_fn, 'rb'))
    if ii == 0:
        for vn in dd.keys():
            if vn == 'ocean_time':
                data_dict[vn] = np.nan * np.ones(NT)
            else:
                data_dict[vn] = np.nan * np.ones(((NT,) + dd[vn].shape))
    for vn in dd.keys():
        ddv = dd[vn]
        if vn == 'ocean_time':
            data_dict[vn][ii] = ddv
        else:
            if ddv.ndim == 2:
                data_dict[vn][ii, :, :] = ddv
            elif ddv.ndim == 3:
                data_dict[vn][ii, :, :, :] = ddv
            else:
                print('problem with the temporary data dict!')
    ii += 1
        
# +++++++++ end parallel subprocess +++++++++++++++++++++++++++++++++

# Write to NetCDF using xarray (this works, hooray!).

tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_clm_file(data_dict, out_fn)
print('\n- Write clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_ini.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_ini_file(in_fn, out_fn)
print('\n- Write ini file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_bry_file(in_fn, out_fn)
print('\n- Write bry file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()
    
def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn, decode_times=False)
    print(ds)
    ds.close()

# check results
nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
if False:
    # print info about the files to the screen
    for fn in nc_list:
        print_info(out_dir / fn)
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
