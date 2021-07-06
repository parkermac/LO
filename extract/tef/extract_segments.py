"""
A tool to extract hourly time series of volume and selected tracers and derived quantities
for budgets in the segments.

To test on mac:
run extract_segments -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -test True

And on perigee:
python extract_segments.py -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -Nproc 20


"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from time import time
import Lfun
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import xarray as xr

import zrfun
import tef_fun

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt00 = time()
print(' Doing segment extraction for '.center(60,'='))
print(' gtagex = ' + Ldir['gtagex'])
outname = 'segments_' + ds0 + '_' + ds1 + '.nc'
print(' outname = ' + outname)
# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
Lfun.make_dir(out_dir)
# name the output file and make sure it does not exist
out_fn = out_dir / outname
out_fn.unlink(missing_ok=True)

# make the scratch directory for holding temporary files
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('segment_temp_' + ds0 + '_' + ds1)
Lfun.make_dir(temp_dir, clean=True)

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
if Ldir['testing']:
    fn_list = fn_list[:3]
else:
    pass
    
print('Doing initial data extraction:')
# We do extractions one hour at a time, as separate subprocess jobs.
# Files are saved to temp_dir.
tt000 = time()
proc_list = []
N = len(fn_list)
for ii in range(N):
    fn = fn_list[ii]
    d = fn.parent.name.replace('f','')
    nhis = int(fn.name.split('.')[0].split('_')[-1])
    cmd_list = ['python3', 'extract_segment_one_time.py',
            '-pth', str(Ldir['roms_out']),
            '-out_dir',str(temp_dir),
            '-gtagex', Ldir['gtagex'],
            '-d', d, '-nhis', str(nhis),
            '-get_bio', str(Ldir['get_bio']),
            '-test', str(Ldir['testing'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    
    Nproc = Ldir['Nproc']
    if (np.mod(ii,Nproc) == 0) or (ii == N-1):
        tt0 = time()
        for proc in proc_list:
            proc.communicate()
        print(' - %d out of %d: %d took %0.2f sec' % (ii, N, Nproc, time()-tt0))
        sys.stdout.flush()
        proc_list = []
print('Total elapsed time = %0.2f sec' % (time()-tt000))

# now load all the results into a single xarray object
A_list = list(temp_dir.glob('A*.p'))
A_list.sort()

# make a list of the datetimes (really should to this in extract_segment_one_time)
dt_list = []
for fn in fn_list:
    ds = nc.Dataset(fn)
    ot = ds['ocean_time'][:]
    dt_list.append(Lfun.modtime_to_datetime(ot.data[0]))

# make a list of the output as DataArrays
x_list = []
for A_fn in A_list:
    A = pickle.load(open(A_fn, 'rb'))
    x_list.append(xr.DataArray(A, dims=('seg','vn')))
# and concatenate into a single DataArray, with time as the concatenating dimension
d = xr.concat(x_list, pd.Index(dt_list, name='time'))
# save it to NetCDF
d.to_netcdf(out_fn)

# Clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()

if Ldir['testing']:
    # check results
    dd = xr.open_dataarray(out_fn)
    ddd = dd.sel(vn='salt', seg='J1')
    D = ddd.values # an ndarray
    print(D)


