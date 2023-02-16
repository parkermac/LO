"""
Process tef2 extractions, giving transport vs. salinity for:
volume, salt, and other variables.

Can be run by making interactive choices in the terminal, or using command line
arguments (as when it is run by extract_sections.py).

PERFORMANCE: ...

To test on mac:
run process_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -test True

"""

import sys
import xarray as xr
import numpy as np
import pickle
from time import time
import pandas as pd

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# import tef_fun

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)


# in_dir00 = Ldir['LOo'] / 'extract'
# if len(args.gtagex) == 0:
#     gtagex = Lfun.choose_item(in_dir00)
# else:
#     gtagex = args.gtagex
# in_dir0 = in_dir00 / gtagex / 'tef'
# if (len(args.ds0)==0) or (len(args.ds1)==0):
#     ext_name = Lfun.choose_item(in_dir0, tag='extractions')
# else:
#     ext_name = 'extractions_' + args.ds0 + '_' + args.ds1
# in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.nc')]
    
# out_dir = in_dir0 / ext_name.replace('extractions', 'processed')
# Lfun.make_dir(out_dir, clean=True)

# vn_list = tef_fun.vn_list # only use extracted variables

# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nProcessing TEF extraction:')
print(str(in_dir))

if Ldir['testing']:
    sect_list = ['mb1.nc']

tt00 = time()

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = xr.open_dataset(in_dir / ext_fn)
    V = dict()
    ds_vn_list = [vn for vn in ds.data_vars]
    for vn in vn_list:
        if vn in ds_vn_list:
            V[vn] = ds[vn].to_numpy()
        else:
            print(' - variable %s not found' % (vn))
    V['salt2'] = V['salt']*V['salt']
    q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy()
    V['q'] = q
    ot = ds['time'].to_numpy()
    zeta = ds['zeta'].to_numpy()
    ds.close()

    # make arrays of property transport
    QV = dict()
    for vn in V.keys():
        if vn == 'q':
            QV[vn] = q
        else:
            QV[vn] = q*V[vn]
    NT, NZ, NX = q.shape
    
    # define salinity bins
    sedges = np.linspace(0, 36, 1001)
    sbins = sedges[:-1] + np.diff(sedges)/2
    NS = len(sbins) # number of salinity bins

    # TEF variables
    Omat = np.zeros((NT, NS))
    TEF = dict()
    for vn in QV.keys():
        TEF[vn] = Omat.copy()
    # other variables
    omat = np.zeros(NT)
    qnet = omat.copy()
    fnet = omat.copy()
    ssh = omat.copy()
    g = 9.8
    rho = 1025

    # process into salinity bins
    for tt in range(NT):
            
        si = V['salt'][tt,:,:].squeeze()
        sf = zfun.fillit(si).flatten()
        sf = sf[~np.isnan(sf)]
        # sort into salinity bins
        inds = np.digitize(sf, sedges, right=True)
        indsf = inds.copy().flatten()
            
        for vn in QV.keys():
            XI = QV[vn][tt,:,:].squeeze()
            XF = zfun.fillit(XI).flatten()
            XF = XF[~np.isnan(XF)]
            if vn == 'q':
                # also keep track of volume transport
                qnet[tt] = XF.sum()
                # and tidal energy flux
                zi = zeta[tt,:].squeeze()
                ssh[tt] = np.nanmean(zfun.fillit(zi))
                fnet[tt] = g * rho * ssh[tt] * qnet[tt]
            counter = 0
            for ii in indsf:
                TEF[vn][tt, ii-1] += XF[counter]
                counter += 1
    
    TEF['ot'] = ot
    TEF['sbins'] = sbins
    TEF['qnet'] = qnet
    TEF['fnet'] = fnet
    TEF['ssh'] = ssh
    pickle.dump(TEF, open(out_dir / out_fn, 'wb'))
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()

print('\nTotal elapsed time = %d seconds' % (time()-tt00))



