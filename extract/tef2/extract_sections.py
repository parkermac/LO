"""
Code to extract tef2 sections.

To test on mac:
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -Nproc 40

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

# vn_list_old = ['salt', 'temp', 'oxygen',
#     'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
#     'TIC', 'alkalinity']
#
# vn_list_new = ['salt', 'temp', 'oxygen',
#     'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN',
#     'TIC', 'alkalinity']

# add custom dict fields
# long_name_dict['q'] = 'transport'
# units_dict['q'] = 'm3 s-1'
# long_name_dict['lon'] = 'longitude'
# units_dict['lon'] = 'degrees'
# long_name_dict['lat'] = 'latitude'
# units_dict['lat'] = 'degrees'
# long_name_dict['h'] = 'depth'
# units_dict['h'] = 'm'
# long_name_dict['z0'] = 'z on rho-grid with zeta=0'
# units_dict['z0'] = 'm'
# long_name_dict['DA0'] = 'cell area on rho-grid with zeta=0'
# units_dict['DA0'] = 'm2'
# long_name_dict['DA'] = 'cell area on rho-grid'
# units_dict['DA'] = 'm2'

# Note: it was not faster to try to parallelize the job this way, but
# maybe it will be more important on the linux machines.

# loop over all jobs
tt0 = time()
N = len(fn_list)
proc_list = []
for ii in range(N):
    # Launch a job and add its process to a list.
    fn = fn_list[ii]
    ii_str = ('0000' + str(ii))[-5:]
    out_fn = temp_dir / ('CC_' + ii_str + '.p')
    cmd_list = ['python3', 'get_one_section.py',
            '-sect_df_fn', str(sect_df_fn),
            '-in_fn',str(fn),
            '-out_fn', str(out_fn),
            '-test', str(Ldir['testing'])]
    
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



