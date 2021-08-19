"""
Process TEF extractions, giving transport vs. salinity for:
volume, salt, salinity-squared and other variables.

Can be run by making interactive choices in the terminal, or using command line
arguments (as when it is run by extract_sections.py).

PERFORMANCE: 24 minutes for a year, all section, all variables, on my mac.

To test on mac:
run process_sections -gtagex cas6_v3_lo8b -0 2019.07.04 -1 2019.07.06

"""

import sys
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import argparse

from lo_tools import Lfun, zfun
import tef_fun
Ldir = Lfun.Lstart()

parser = argparse.ArgumentParser()
parser.add_argument('-gtagex', type=str, default='')   # e.g. cas6_v3_lo8b
parser.add_argument('-0', '--ds0', type=str, default='')        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='')        # e.g. 2019.07.04
args = parser.parse_args()

in_dir00 = Ldir['LOo'] / 'extract'
if len(args.gtagex) == 0:
    gtagex = Lfun.choose_item(in_dir00)
else:
    gtagex = args.gtagex
in_dir0 = in_dir00 / gtagex / 'tef'
if (len(args.ds0)==0) or (len(args.ds1)==0):
    ext_name = Lfun.choose_item(in_dir0, tag='extractions')
else:
    ext_name = 'extractions_' + args.ds0 + '_' + args.ds1
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.nc')]
    
out_dir = in_dir0 / ext_name.replace('extractions', 'processed')
Lfun.make_dir(out_dir, clean=True)

vn_list = tef_fun.vn_list # only use extracted variables

print('\nProcessing TEF extraction:')
print(str(in_dir))

tt00 = time()

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = nc.Dataset(in_dir / ext_fn)
    V = dict()
    ds_vn_list = [vn for vn in ds.variables]
    for vn in vn_list:
        if vn in ds_vn_list:
            V[vn] = ds[vn][:]
        else:
            print(' - variable %s not found' % (vn))
    V['salt2'] = V['salt']*V['salt']
    q = ds['q'][:]
    V['q'] = q
    ot = ds['ocean_time'][:]
    zeta = ds['zeta'][:]
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



