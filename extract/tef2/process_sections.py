"""
Process tef2 extractions, giving transport vs. salinity for:
volume, salt, and other variables.

Can be run by making interactive choices in the terminal, or using command line
arguments (as when it is run by extract_sections.py).

PERFORMANCE: Using the new "binned_statistic" method this is significantly faster,
about 5 minutes for 83 section for a year with only salt. Is this really faster?

To test on mac:
run process_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True
run process_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06

And for a full year:

(this only has salt)
run process_sections -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True

"""

import sys
import xarray as xr
import numpy as np
import pickle
from time import time
import pandas as pd
from scipy.stats import binned_statistic

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

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    sect_list = ['jdf3.nc']
    
# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nProcessing TEF extraction:')
print(str(in_dir))

tt00 = time()

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = xr.open_dataset(in_dir / ext_fn)
    V = dict()
    for vn in vn_list:
        V[vn] = ds[vn].to_numpy()
    # V['salt2'] = V['salt']*V['salt']
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
    if Ldir['testing']:
        NS = 36 # number of salinity bins
    else:
        NS = 1000 # number of salinity bins
    S_low = 0
    S_hi = 36
    sedges = np.linspace(S_low, S_hi, NS+1)
    sbins = sedges[:-1] + np.diff(sedges)/2

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

    # Process into salinity bins.
    # NOTE: this seems like a lot of nested loops (section/time/variable/digitization indices)
    # but it runs reasonably fast. Also we only have to run it once.
    for tt in range(NT):
            
        if False:
            sf = V['salt'][tt,:,:].squeeze().flatten()
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
                    ssh[tt] = zi.mean()
                    fnet[tt] = g * rho * ssh[tt] * qnet[tt]
                counter = 0
                for ii in indsf:
                    TEF[vn][tt, ii-1] += XF[counter]
                    counter += 1
                
                
                if Ldir['testing'] and (tt==10) and (vn=='salt'):
                    print(TEF[vn][tt,:])
                
        else:
            # alternate version
            
            sf = V['salt'][tt,:,:].squeeze().flatten()
            
            for vn in QV.keys():
                XF = QV[vn][tt,:,:].squeeze().flatten()
                if vn == 'q':
                    # also keep track of volume transport
                    qnet[tt] = XF.sum()
                    # and tidal energy flux
                    zi = zeta[tt,:].squeeze()
                    ssh[tt] = zi.mean()
                    fnet[tt] = g * rho * ssh[tt] * qnet[tt]
                    
                # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
                TEF[vn][tt,:] = binned_statistic(sf, XF, statistic='sum', bins=NS, range=(S_low,S_hi)).statistic
                
                # results are identical
                
                
                if Ldir['testing'] and (tt==10) and (vn=='salt'):
                    print(TEF[vn][tt,:])
            
        
    
    TEF['ot'] = ot
    TEF['sbins'] = sbins
    TEF['qnet'] = qnet
    TEF['fnet'] = fnet
    TEF['ssh'] = ssh
    pickle.dump(TEF, open(out_dir / out_fn, 'wb'))
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()

print('\nTotal elapsed time = %d seconds' % (time()-tt00))



