"""
Code to extract tef2 sections.

To test on mac:
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -test True
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -Nproc 10 -get_bio True

Doing this with subprocesses (Nproc = 10, on my mac) was about 2x as fast as doing
it sequentially within this program. The less-than-expected speedup may be because
each job only takes about a second, and there is overhead to spinning up new python
jobs because of imports.

Also, this is a memory-intensive calculation, so be careful about using Nproc > 10
(10 is the default in extract_argfun).

Performance: took about 1-2 sec per history file (Nproc = 10, on my mac).
- 58 sec per day with get_bio True (11 3-D variables)
- 24 sec per day with get_bio False (only salt)
- a bit over an hour per year on apogee (only salt)

"""

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import sys
import pandas as pd
import xarray as xr
import numpy as np
import pickle

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

fn_list = Lfun.get_fn_list('hourly', Ldir, Ldir['ds0'], Ldir['ds1'])

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
out_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
temp_dir = out_dir0 / ('temp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_dir, clean=True)

if Ldir['testing']:
    fn_list = fn_list[:3]
if Ldir['get_bio']:
    vn_type = 'bio'
else:
    vn_type = 'salt'

# loop over all jobs
tt0 = time()
N = len(fn_list)
proc_list = []
for ii in range(N):
    # Launch a job and add its process to a list.
    fn = fn_list[ii]
    ii_str = ('0000' + str(ii))[-5:]
    out_fn = temp_dir / ('CC_' + ii_str + '.nc')
    # use subprocesses
    cmd_list = ['python3', 'get_one_section.py',
            '-sect_df_fn', str(sect_df_fn),
            '-in_fn',str(fn),
            '-out_fn', str(out_fn),
            '-vn_type', vn_type]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    # If we have accumulated Nproc jobs, or are at the end of the
    # total number of jobs, then stop and make sure all the jobs
    # in proc_list have finished, using the communicate method.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                print('\nSTDOUT:')
                print(stdout.decode())
                sys.stdout.flush()
            if len(stderr) > 0:
                print('\nSTDERR:')
                print(stderr.decode())
                sys.stdout.flush()
        # Then initialize a new list.
        proc_list = []
    # Print screen output about progress.
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == N-1):
        print(str(ii))
        sys.stdout.flush()
    
print('Total processing time = %0.2f sec' % (time()-tt0))

# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','CC'], stdin=pp1.stdout, stdout=Pi)
temp_fn = str(temp_dir)+'/all.nc'
cmd_list = ['ncrcat','-p', str(temp_dir), '-O', temp_fn]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\nSTDOUT:')
    print(stdout.decode())
    sys.stdout.flush()
if len(stderr) > 0:
    print('\nSTDERR:')
    print(stderr.decode())
    sys.stdout.flush()
        
"""
Next we want to repackage these results into one NetCDF file per section, with all times.

We will mostly follow the structure of the output of LO/tef/extract_sections.py
so that we can mostly recycle the subsequent processing code:

The result, looking for example at "this_ds" created below
for one section while testing:
    
<xarray.Dataset>
Dimensions:  (time: 3, p: 8, z: 30)
Coordinates:
  * time     (time) datetime64[ns] 2021-07-04 ... 2021-07-04T02:00:00
Dimensions without coordinates: p, z
Data variables:
    h        (p) float64 17.14 20.32 22.74 21.58 20.65 18.38 16.7 14.99
    dd       (p) float64 ...
    zeta     (time, p) float32 0.2535 0.2531 0.2528 ... 0.0301 0.03066 0.03053
    salt     (time, z, p) float32 ...
    vel      (time, z, p) float64 ...
    DZ       (time, z, p) float64 0.7218 0.864 0.9729 ... 0.2256 0.217 0.2074

The dimension "p" means a point on the stairstep section. "dd" is the point width [m],
and "DZ" [m] is the vertical thickness of each cell.
    
"""

ds1 = xr.open_dataset(temp_fn)
S = zrfun.get_basic_info(fn_list[0], only_S=True)
eta = ds1.zeta.values.squeeze() # packed (t, p)
NT, NP = eta.shape
hh = ds1.h.values.squeeze().reshape(1,NP) * np.ones((NT,1))
zw = zrfun.get_z(hh, eta, S, only_w=True)
dz = np.diff(zw, axis=0) # NOTE: this is packed (z,t,p)
DZ = np.transpose(dz, (1,0,2)) # packed (t,z,p)

sect_list = list(sect_df.sn.unique())
sect_list.sort()
for sn in sect_list:
    """
    A useful tool for pulling out a section is np.where() combined with the
    xr Dataset method isel(), as is done here.
    """
    ii = np.where(sect_df.sn == sn)[0]
    this_ds = ds1.isel(p=ii)
    # add DZ
    this_DZ = DZ[:,:,ii]
    this_ds['DZ'] = (('time','z','p'), this_DZ)
    this_fn = out_dir / (sn + '.nc')
    this_ds.to_netcdf(this_fn)
    
ds1.close()


