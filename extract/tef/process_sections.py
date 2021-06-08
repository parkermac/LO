"""
Process TEF extractions, giving transport vs. salinity for:
volume, salt, and salinity-squared.

Also velocity and area which are relevant to some tidal-pumping calculations
I am trying (that may not lead anywhere).

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

sect_list = [item.name for item in in_dir.glob('*') if '.nc' in item.name]
    
out_dir = in_dir0 / ext_name.replace('extractions', 'processed')
Lfun.make_dir(out_dir, clean=True)

for ext_fn in sect_list:
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = nc.Dataset(in_dir / ext_fn)
    q = ds['q'][:]
    da = ds['DA'][:] # new
    s = ds['salt'][:]
    ot = ds['ocean_time'][:]
    zeta = ds['zeta'][:]
    ds.close()

    # TEF sort into salinity bins
    qs = q*s
    qs2 = q*s*s
    NT, NZ, NX = q.shape
    # initialize intermediate results arrays for TEF quantities
    sedges = np.linspace(0, 36, 1001) # original was 1001 used 5001 for Willapa
    sbins = sedges[:-1] + np.diff(sedges)/2
    NS = len(sbins) # number of salinity bins

    # TEF variables
    tef_q = np.zeros((NT, NS))
    tef_vel = np.zeros((NT, NS))
    tef_da = np.zeros((NT, NS))
    tef_qs = np.zeros((NT, NS))
    tef_qs2 = np.zeros((NT, NS))

    # other variables
    qnet = np.zeros(NT)
    fnet = np.zeros(NT)
    ssh = np.zeros(NT)
    g = 9.8
    rho = 1025

    for tt in range(NT):
        if np.mod(tt,1000) == 0:
            print('  time %d out of %d' % (tt,NT))
            sys.stdout.flush()
            
        qi = q[tt,:,:].squeeze()
        if isinstance(qi, np.ma.MaskedArray):
            qf = qi[qi.mask==False].data.flatten()
        else:
            qf = qi.flatten()
            
        si = s[tt,:,:].squeeze()
        if isinstance(si, np.ma.MaskedArray):
            sf = si[qi.mask==False].data.flatten()
        else:
            sf = si.flatten()
            
        dai = da[tt,:,:].squeeze()
        if isinstance(dai, np.ma.MaskedArray):
            daf = dai[qi.mask==False].data.flatten()
        else:
            daf = dai.flatten()
            
        qsi = qs[tt,:,:].squeeze()
        if isinstance(qsi, np.ma.MaskedArray):
            qsf = qsi[qi.mask==False].data.flatten()
        else:
            qsf = qsi.flatten()
            
        qs2i = qs2[tt,:,:].squeeze()
        if isinstance(qs2i, np.ma.MaskedArray):
            qs2f = qs2i[qi.mask==False].data.flatten()
        else:
            qs2f = qs2i.flatten()
            
        # sort into salinity bins
        inds = np.digitize(sf, sedges, right=True)
        indsf = inds.copy().flatten()
        counter = 0
        for ii in indsf:
            tef_q[tt, ii-1] += qf[counter]
            tef_da[tt, ii-1] += daf[counter] # new
            tef_qs[tt, ii-1] += qsf[counter]
            tef_qs2[tt, ii-1] += qs2f[counter]
            counter += 1
        
        # also keep track of volume transport
        qnet[tt] = qf.sum()
        
        # and tidal energy flux
        zi = zeta[tt,:].squeeze()
        ff = zi.reshape((1,NX)) * qi
        fnet[tt] = g * rho * ff.sum()

    # save results
    tef_dict = dict()
    tef_dict['tef_q'] = tef_q
    tef_dict['tef_da'] = tef_da
    tef_dict['tef_qs'] = tef_qs
    tef_dict['tef_qs2'] = tef_qs2
    tef_dict['sbins'] = sbins
    tef_dict['ot'] = ot
    tef_dict['qnet'] = qnet
    tef_dict['fnet'] = fnet
    tef_dict['ssh'] = np.mean(zeta, axis=1)
    pickle.dump(tef_dict, open(out_dir / out_fn, 'wb'))
    



