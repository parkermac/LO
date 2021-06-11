"""
Process TEF extractions, giving transport vs. salinity for:
volume, salt, and salinity-squared.

Takes about 30 minutes per year for all 39 cas6 sections, so you might want to
do ctrl-z, bg, to send it to background after making the choices.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle

import Lfun
import zfun

Ldir = Lfun.Lstart()

in_dir00 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir00)
in_dir0 = in_dir00 / gtagex / 'tef'
ext_name = Lfun.choose_item(in_dir0, tag='extractions')
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.nc')]
    
out_dir = in_dir0 / ext_name.replace('extractions', 'processed')
Lfun.make_dir(out_dir, clean=True)

vn_list = ['salt', 'temp', 'oxygen', 'NO3', 'TIC', 'alkalinity', 'q']

for ext_fn in sect_list:
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = nc.Dataset(in_dir / ext_fn)
    
    
    
    V = dict()
    for vn in vn_list:
        if vn in ds.variables:
            V[vn] = ds[vn][:]
        else:
            print('variable %s not found' % (vn))
            vn_list.remove(vn)
    V['salt2'] = V['salt']*V['salt']
    
    q = ds['q'][:]
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

    for tt in range(NT):
        if np.mod(tt,1000) == 0:
            print('  time %d out of %d' % (tt,NT))
            sys.stdout.flush()
            
        si = V['salt'][tt,:,:].squeeze()
        if isinstance(si, np.ma.MaskedArray):
            sf = si[si.mask==False].data.flatten()
        else:
            sf = si.flatten()
        sf = sf[~np.isnan(sf)]
        # sort into salinity bins
        inds = np.digitize(sf, sedges, right=True)
        indsf = inds.copy().flatten()
            
        for vn in QV.keys():
            XI = QV[vn][tt,:,:].squeeze()
            if isinstance(XI, np.ma.MaskedArray):
                XF = XI[XI.mask==False].data.flatten()
            else:
                XF = XI.flatten()
            XF = XF[~np.isnan(XF)]
                
            if vn == 'q':
                # also keep track of volume transport
                qnet[tt] = XF.sum()
                # and tidal energy flux
                zi = zeta[tt,:].squeeze()
                ff = zi.reshape((1,NX)) * XI
                fnet[tt] = g * rho * ff.sum()
                
            
            counter = 0
            for ii in indsf:
                TEF[vn][tt, ii-1] += XF[counter]
                counter += 1
    
    TEF['ot'] = ot
    TEF['sbins'] = sbins
    TEF['qnet'] = qnet
    TEF['fnet'] = fnet
    TEF['ssh'] = np.mean(zeta, axis=1)
    
    pickle.dump(TEF, open(out_dir / out_fn, 'wb'))
    



