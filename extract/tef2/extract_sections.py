"""
Code to extract tef2 sections.

To test on mac:
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -Nproc 10 -get_bio True

Doing this with subprocesses (Nproc = 10, on my mac) was about 2x as fast as doing
it sequentially within this program. The less-than-expected speedup may be because
each job only takes about a second, and there is overhead to spinning up new python
jobs because of imports.

Also, this is a memory-intensive calculation, so be careful about using Nproc > 10
(10 is the default in extract_argfun).

Performance: took about 1-2 sec per history file (Nproc = 10, on my mac).
- 58 sec per day with get_bio True (11 3-D variables)
- 24 sec per day with get_bio False (only salt)

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
    fn_list = [fn_list[0]]
    
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
    out_fn = temp_dir / ('CC_' + ii_str + '.p')
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

"""
Next we want to repackage these results into one NetCDF file per section, with all times.

We will follow the structure of the output of LO/tef/extract_sections.py so that we can mostly
recycle the subsequent processing code:

Variables in the NetCDF files:
- salt is hourly salinity in each cell (t, z, x-or-y) [same for all other variables]
- q is hourly transport in each cell (t, z, x-or-y)
- vel is velocity in each cell (t, z, x-or-y) positive to East or North
- DA is the area of each cell (t, z, x-or-y) hence: q = vel * DA
- z0 is the average z-position of cell centers (assumes SSH=0), useful for plotting
- DA0 is the average cross-sectional area of each cell (assumes SSH=0)
- h is depth on the section (x-or-y) positive down
- zeta is SSH on the section (t, x-or-y) positive up
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.
"""

sect_list = list(sect_df.sn.unique())
sect_list.sort()
cc_list = list(temp_dir.glob('CC_*.p'))
cc_list.sort()

S = zrfun.get_basic_info(fn_list[0], only_S=True)
NZ = S['N']
NT = len(cc_list) # number of times
for sn in ['mb9']: sect_list:
    df = sect_df[sect_ds.sn==sn]
    ii = df.index.to_numpy()
    NX = len(ii)
    
    # initialize arrays
    
    # initialize DataSet
    
    # loop over all times to fill arrays
    
    # then add these to the DataSet and save to output.

