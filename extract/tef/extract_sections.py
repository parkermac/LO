"""
Extract fields at a number of sections which may be used later for TEF analysis
of transport and transport-weighted properties, making use of multiple subprocesses
to speed up operation.

All input parameters are specified at the command line, so this can be run in the background.

Performance: 4 hours per year, perigee, 39 cas6 sections, -get_bio True.  FAST!

To test on mac (default is to just get ai1 section):

run extract_sections -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -test True

To get all sections for the same time:

python extract_sections -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.04 -get_bio True -sect_name all > log_all_test &

python extract_sections.py -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -sect_name all > log_all_test &

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
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import zrfun
import tef_fun
if Ldir['testing']:
    from importlib import reload
    reload(tef_fun)
    
# set list of variables to extract
if Ldir['get_bio']:
    vn_list = tef_fun.vn_list
else:
    vn_list = ['salt']

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']
dt0 = datetime.strptime(ds0, Ldir['ds_fmt'])
dt1 = datetime.strptime(ds1, Ldir['ds_fmt'])
ndays = (dt1-dt0).days + 1

tt00 = time()
print('Working on:')
outname = 'extractions_' + ds0 + '_' + ds1
print(outname)

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / outname
Lfun.make_dir(out_dir, clean=True)

# make the scratch directory
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('tef_temp_' + ds0 + '_' + ds1)
Lfun.make_dir(temp_dir, clean=True)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
# initialize a dictionary of info for each section
sect_info = dict()
# select which sections to extract
if Ldir['sect_name'] == 'all':
    # full list
    sect_list = [item for item in sect_df.index]
else: # single item
    if Ldir['sect_name'] in sect_df.index:
        sect_list = [Ldir['sect_name']]
    else:
        print('That section is not available')
        sys.exit()

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
NT = len(fn_list)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
NZ = S['N']

# Create and save the sect_info dict
tt0 = time()
# - make a dictionary of info for each section
sect_info = dict()
print('\nGetting section definitions and indices:')
for sect_name in sect_list:
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    # - get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    NX = len(Mask)
    # - save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, NX, Lon, Lat)
info_fn = temp_dir / 'sect_info.p'
pickle.dump(sect_info, open(info_fn, 'wb'))
print('Elapsed time = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
print('\nInitializing NetCDf output files:')
sl = sect_list.copy()
while len(sl) > 0:
    if len(sl) >= 10:
        print(sl[:10])
        sl = sl[10:]
    else:
        print(sl)
        sl = []
for sect_name in sect_list:
    out_fn = out_dir / (sect_name + '.nc')
    ii0, ii1, jj0, jj1, sdir, NX, Lon, Lat = sect_info[sect_name]
    tef_fun.start_netcdf(fn, out_fn, NT, NX, NZ, Lon, Lat, Ldir, vn_list)
print('Elapsed time = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# do the initial data extraction
# if Ldir['testing']:
#     fn_list = [fn_list[0]]

print('\nDoing initial data extraction:')
# We do extractions one hour at a time, as separate subprocess jobs.
# Running Nproc (e.g. 20) of these in parallel makes the code much faster.
tt000 = time()
proc_list = []
N = len(fn_list)
for ii in range(N):
    fn = fn_list[ii]
    d = fn.parent.name.replace('f','')
    nhis = int(fn.name.split('.')[0].split('_')[-1])
    cmd_list = ['python3', 'extract_one_time.py',
            '-pth', str(Ldir['roms_out']),
            '-out_dir',str(temp_dir),
            '-gtagex', Ldir['gtagex'],
            '-d', d, '-nhis', str(nhis),
            '-get_bio', str(Ldir['get_bio'])]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    
    Nproc = 20
    if (np.mod(ii,Nproc) == 0) or (ii == N-1):
        tt0 = time()
        for proc in proc_list:
            proc.communicate()
        print(' - %d out of %d: %d took %0.2f sec' % (ii, N, Nproc, time()-tt0))
        sys.stdout.flush()
        proc_list = []
print('Elapsed time = %0.2f sec' % (time()-tt000))

# write fields to NetCDF

# # extract and save time-dependent fields
tt0 = time()
print('\nAdding data to NetCDF files:')
for sect_name in sect_list:
    out_fn = out_dir / (sect_name + '.nc')
    # this is where we add the data from this history file
    # to all of the sections, each defined by sinfo
    tef_fun.add_fields(out_fn, temp_dir, sect_name, vn_list, S, NT)
print('Elapsed time = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
    
print('\nTotal elapsed time = %d seconds' % (time()-tt00))

# clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()
    


