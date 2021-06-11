"""
This code extracts all needed TEF data for one history file,
looping over all variables and all sections.
"""
from pathlib import Path
import sys
from datetime import datetime, timedelta
import argparse
import numpy as np
import netCDF4 as nc
from time import time

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import zfun, zrfun

import tef_fun

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pth', type=str, default='/Users/pm8/Documents/LiveOcean_roms/output')
parser.add_argument('-gtagex', type=str, default='cas6_v3_lo8b')
parser.add_argument('-d', type=str, default='2019.07.04')
parser.add_argument('-nhis', type=int, default=1)
parser.add_argument('-get_bio', type=Lfun.boolean_string, default=True)
args = parser.parse_args()

nhiss = ('0000' + str(args.nhis))[-4:]
fn = Path(args.pth) / args.gtagex / ('f' + args.d) / ('ocean_his_' + nhiss + '.nc')

# set list of variables to extract
if args.get_bio:
    vn_list = ['salt', 'temp', 'oxygen', 'NO3', 'TIC', 'alkalinity']
else:
    vn_list = ['salt']

tt0 = time()
# Make the sect_info dict (could be done by calling function and then just loaded here)
# - get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
sect_list = [item for item in sect_df.index]
# - get grid info
G = zrfun.get_basic_info(fn, only_G=True)
# - make a dictionary of info for each section
sect_info = dict()
print('\nGetting section definitions and indices:')
for sect_name in sect_list:
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    # - get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    # - save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir)
print('Time to get sect_info = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# gather fields for section extractions
tt0 = time()
ds = nc.Dataset(fn)
Zeta = ds['zeta'][0,:,:]
U = ds['u'][0,:,:,:]
V = ds['v'][0,:,:,:]
CC = dict()
for vn in vn_list:
    CC[vn] = ds[vn][0,:,:,:]
ds.close()
print('Time to gather fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# do section extractions
tt0 = time()
A = dict()
for sect_name in sect_list:
    ii0, ii1, jj0, jj1, sdir = sect_info[sect_name]
    if sdir=='NS':
        h = G['h'][jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=1)
        zeta = Zeta[jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=1)
        d = G['DY'][jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=1)
    elif sdir=='EW':
        h = G['h'][jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=0)
        zeta = Zeta[jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=0)
        d = G['DX'][jj0:jj1+1,ii0:ii1+1].squeeze().mean(axis=0)
    # then velocity
    if sdir=='NS':
        vel = U[:, jj0:jj1+1, ii0].squeeze()
    elif sdir=='EW':
        vel = V[:, jj0, ii0:ii1+1].squeeze()
    # save the tracer fields averaged onto this section
    C = dict()
    for vn in vn_list:
        if sdir=='NS':
            c = (CC[vn][:,jj0:jj1+1,ii0].squeeze()
                + CC[vn][:,jj0:jj1+1,ii1].squeeze())/2
        elif sdir=='EW':
            c = (CC[vn][:,jj0,ii0:ii1+1].squeeze()
                + CC[vn][:,jj1,ii0:ii1+1].squeeze())/2
        C[vn] = zfun.fillit(c)
    C['h'] = zfun.fillit(h)
    C['zeta'] = zfun.fillit(zeta)
    C['d'] = zfun.fillit(d)
    C['vel'] = zfun.fillit(vel)
    
    A[sect_name] = C
print('Time to extract all section data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

