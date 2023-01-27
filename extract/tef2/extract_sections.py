"""
Code to extract tef2 sections.

To test on mac:
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True

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
sect_df = pd.read_pickle(tef2_dir / ('sect_df_' + gctag + '.p'))

fn_list = Lfun.get_fn_list('hourly', Ldir, Ldir['ds0'], Ldir['ds1'])

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
out_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
temp_dir = out_dir0 / ('temp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_dir, clean=True)

if Ldir['testing']:
    sect_df = sect_df.loc[(sect_df.sn == 'mb8') | (sect_df.sn == 'mb9'),:].copy()
    sect_df = sect_df.reset_index(drop=True)
    
vn_list = ['salt']

if True:#Ldir['testing']:
    fn_list = [fn_list[0]]

# grid info
fn = fn_list[0]
S = zrfun.get_basic_info(fn, only_S=True)
ds = xr.open_dataset(fn)
DX = 1/ds.pm.values
DY = 1/ds.pn.values
# Get spacing on u and v grids
# dxu = DX[:,:-1] + np.diff(DX,axis=1)/2
# dyv = DY[:-1,:] + np.diff(DY,axis=0)/2
dxv = DX[:-1,:] + np.diff(DX,axis=0)/2
dyu = DY[:,:-1] + np.diff(DY,axis=1)/2
    
out_df = pd.DataFrame()

u_df = sect_df[sect_df.uv == 'u']
v_df = sect_df[sect_df.uv == 'v']

# Fields that do not change with time
C = dict()
CC = dict()
h = ds.h.values
CC['h'] = (h[sect_df.jrp, sect_df.irp]  + h[sect_df.jrm, sect_df.irm])/2
dxvv = dxv[v_df.j, v_df.i]
dyuu = dyu[u_df.j, u_df.i]
dd = np.nan * np.ones(CC['h'].shape)
dd[v_df.index] = dxvv
dd[u_df.index] = dyuu

tt0 = time()
for fn in fn_list:
    # eventually this will be handled as simultaneous subprocess jobs
    ds = xr.open_dataset(fn)
    print(fn.name)
    
    # First: tracers and zeta
    for vn in vn_list:
        C[vn] = ds[vn].values.squeeze()
    C['zeta'] = ds.zeta.values.squeeze()
    for vn in vn_list:
        CC[vn] = (C[vn][:, sect_df.jrp, sect_df.irp]  + C[vn][:, sect_df.jrm, sect_df.irm])/2
    CC['zeta'] = (C['zeta'][sect_df.jrp, sect_df.irp]  + C['zeta'][sect_df.jrm, sect_df.irm])/2
    
    # Then: velocity
    u = ds.u.values.squeeze()
    v = ds.v.values.squeeze()
    uu = u[:, u_df.j, u_df.i] * u_df.pm.to_numpy().reshape(1,-1)
    vv = v[:, v_df.j, v_df.i] * v_df.pm.to_numpy().reshape(1,-1)
    # then merge these back into one
    vel = np.nan * np.ones(CC['salt'].shape)
    # I love fancy indexing!
    vel[:,u_df.index] = uu
    vel[:,v_df.index] = vv
    
print('elapsed time = %0.1f sec' % (time()-tt0))

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

sys.exit()

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



