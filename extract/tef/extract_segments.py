"""
A tool to extract hourly time series of volume and selected tracers and derived quantities
for budgets in the segments.

To test on mac:
run extract_segments -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -test True

And on perigee:
python extract_segments.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -Nproc 20 > test2019.log &

Performance: about 4 hours per year on perigee with bio.
"""
from lo_tools import Lfun, zrfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing
import tef_fun

import sys
from time import time
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import xarray as xr

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
    if ((np.mod(ii,Nproc) == 0) and (ii > 0)) or (ii == N-1):
        tt0 = time()
        for proc in proc_list:
            proc.communicate()
        print(' - %d out of %d: %d took %0.2f sec' % (ii, N, Nproc, time()-tt0))
        sys.stdout.flush()
        proc_list = []
print('Total elapsed time = %0.2f sec' % (time()-tt000))

"""
Now we repackage all these single-time extractions into an xarray Dataset.  This is
a little tricky because for the Dataset we want to have it composed of DataArrays that
are variable(time, segment), but what we are starting from is a collection of
time(segment, variable) pandas DataFrames.

We will proceed by using a two step process, first concatenating all the pandas DataFrames
into an xarray DataArray with dimensions(time, segment, variable).  Then we will
convert this to an xarray Dataset with a collection of variable(time,segment) DataArrays.

I'm sure there is a more clever way to do this in xarray, but I am not yet proficient
enough with that module.
"""

# get a list of all our pandas DataFrames
A_list = list(temp_dir.glob('A*.p'))
A_list.sort()
# make a list of the datetimes and form a time index
ot_list = []
for fn in fn_list:
    ds = xr.open_dataset(fn)
    ot = ds['ocean_time'].values[0]
    ot_list.append(ot)
ot_ind = pd.Index(ot_list)
# make a list of the output as DataArrays
x_list = []
for A_fn in A_list:
    A = pd.read_pickle(A_fn)
    x_list.append(xr.DataArray(A, dims=('seg','vn')))
# and concatenate that list into a single DataArray, with time as the concatenating dimension
da = xr.concat(x_list, pd.Index(ot_ind, name='time'))
# repackage the DataArray as a Dataset
vns = da.coords['vn'].values
segs = da.coords['seg'].values
times = da.coords['time'].values
ds = xr.Dataset(coords={'time': times,'seg': segs})
for vn in vns:
    v = da.sel(vn=vn).values
    ds[vn] = (('time','seg'), v)
# save it to NetCDF
ds.to_netcdf(out_fn)
ds.close()

# Clean up
if not Ldir['testing']:
    Lfun.make_dir(temp_dir, clean=True)
    temp_dir.rmdir()

if Ldir['testing']:
    # check results
    dd = xr.open_dataset(out_fn)
    print(dd.salt.sel(seg='J1').values)
    dd.close()


