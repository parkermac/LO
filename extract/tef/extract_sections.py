"""
Extract fields at a number of sections which may be used later for TEF analysis
of transport and other properties, making use of multiple subprocesses
to speed up operation.  It also runs the process_sections.py and bulk_calc.py jobs
unless you use -test True.  You need to run for at least three days to get any
results at the end, because of the tidal averaging.

All input parameters are specified at the command line, so this can be run in the background.
This is essential for year-long extractions.

PERFORMANCE: 6 hours per year, perigee, cas6_v3_lo8b, 39 sections, all variables.  About 6.5 hours
including the processing and bulk_calc steps.  2.25 hours when just extracting salt.

To test on mac (default is to just get salt on section ai1):
run extract_sections -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06

To get all sections and all variables use these flags:
-get_bio True -sect_name all

On perigee use -Nproc 20

"""

from lo_tools import Lfun, zrfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing
import tef_fun
if Ldir['testing']:
    from importlib import reload
    reload(tef_fun)

import sys
from time import time
import numpy as np
import netCDF4 as nc
import pickle
from subprocess import Popen as Po
from subprocess import PIPE as Pi
    
# set list of variables to extract
if Ldir['get_bio']:
    vn_list = tef_fun.vn_list
else:
    vn_list = ['salt']

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt00 = time()
print(' Doing TEF extraction for '.center(60,'='))
print(' gtagex = ' + Ldir['gtagex'])
outname = 'extractions_' + ds0 + '_' + ds1
print(' outname = ' + outname)

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / outname
Lfun.make_dir(out_dir, clean=True)

# make the scratch directory for holding temporary files
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('tef_temp_' + ds0 + '_' + ds1)
Lfun.make_dir(temp_dir, clean=True)

# get the DataFrame of all sections
gridname=Ldir['gtagex'].split('_')[0]
sect_df = tef_fun.get_sect_df(gridname)

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
# - make a dictionary of info for each section
sect_info = dict()
for sect_name in sect_list:
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    # - get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    NX = len(Mask)
    # - save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, NX, Lon, Lat)
info_fn = temp_dir / 'sect_info.p'
pickle.dump(sect_info, open(info_fn, 'wb'))

# Initialize NetCDF output files
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

print('Doing initial data extraction:')
# We do extractions one hour at a time, as separate subprocess jobs.
# Running Nproc (e.g. 20) of these in parallel makes the code much faster.
# Files are saved to temp_dir.
tt000 = time()
proc_list = []
N = len(fn_list)
for ii in range(N):
    fn = fn_list[ii]
    d = fn.parent.name.replace('f','')
    nhis = int(fn.name.split('.')[0].split('_')[-1])
    cmd_list = ['python3', 'extract_section_one_time.py',
            '-pth', str(Ldir['roms_out']),
            '-out_dir',str(temp_dir),
            '-gtagex', Ldir['gtagex'],
            '-d', d, '-nhis', str(nhis),
            '-get_bio', str(Ldir['get_bio'])]
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
print('Elapsed time = %0.2f sec' % (time()-tt000))

# Extract and save time-dependent fields
for sect_name in sect_list:
    out_fn = out_dir / (sect_name + '.nc')
    tef_fun.add_fields(out_fn, temp_dir, sect_name, vn_list, S, NT)
    
# Clean up
Lfun.make_dir(temp_dir, clean=True)
temp_dir.rmdir()

# Then do the processing and bulk calculation
if Ldir['testing'] == False:
    
    # processing
    tt0 = time()
    print(' process_sections '.center(60,'='))
    cmd_list = ['python3', 'process_sections.py',
            '-gtagex', Ldir['gtagex'],
            '-0', ds0, '-1', ds1]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    #print(stdout.decode())
    if len(stderr) > 0:
        print(stderr.decode())
    print('Elapsed time = %0.2f sec' % (time()-tt0))
    
    # bulk_calc
    tt0 = time()
    print(' bulk_calc '.center(60,'='))
    cmd_list = ['python3', 'bulk_calc.py',
            '-gtagex', Ldir['gtagex'],
            '-0', ds0, '-1', ds1]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    #print(stdout.decode())
    if len(stderr) > 0:
        print(stderr.decode())
    print('Elapsed time = %0.2f sec' % (time()-tt0))

print(' Total elapsed time = %d sec '.center(60,'-') % (time()-tt00))
print(' DONE '.center(60,'='))

    


