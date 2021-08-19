"""
This code extracts all needed TEF data for one history file,
looping over all variables and all sections.

Performance: 2-3 sec on perigee.
"""
from pathlib import Path
import sys
import argparse
import numpy as np
import netCDF4 as nc
from time import time
import pickle

from lo_tools import Lfun, zfun, zrfun
import tef_fun

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pth', type=str)
parser.add_argument('-out_dir', type=str)
parser.add_argument('-gtagex', type=str)
parser.add_argument('-d', type=str) # date string like 2019.07.04
parser.add_argument('-nhis', type=int, default=1) # history file number, 1-25
parser.add_argument('-get_bio', type=Lfun.boolean_string, default=True)
parser.add_argument('-testing', type=Lfun.boolean_string, default=False)
args = parser.parse_args()

nhiss = ('0000' + str(args.nhis))[-4:]
fn = Path(args.pth) / args.gtagex / ('f' + args.d) / ('ocean_his_' + nhiss + '.nc')
out_dir = Path(args.out_dir)
Lfun.make_dir(out_dir)
out_fn  = out_dir / ('A_' + args.d + '_' + nhiss + '.p')
out_fn.unlink(missing_ok=True)

G = zrfun.get_basic_info(fn, only_G=True)

# set list of variables to extract
if args.get_bio:
    vn_list = tef_fun.vn_list
else:
    vn_list = ['salt']

if args.testing == True:
    # Make the sect_info dict if we ask it to (otherwise made by calling function)
    tt0 = time()
    # - get the DataFrame of all sections
    gridname=args.gtagex.split('_')[0]
    sect_df = tef_fun.get_sect_df(gridname)
    sect_list = [item for item in sect_df.index]
    # - get grid info
    # - make a dictionary of info for each section
    sect_info = dict()
    print('\nGetting section definitions and indices:')
    for sect_name in sect_list:
        x0, x1, y0, y1 = sect_df.loc[sect_name,:]
        # - get indices for this section
        ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
        NX = len(Mask)
        # - save some things for later use
        sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, NX, Lon, Lat)
    print('Time to get sect_info = %0.2f sec' % (time()-tt0))
    sys.stdout.flush()
else:
    info_fn = out_dir / 'sect_info.p'
    if info_fn.is_file():
        sect_info = pickle.load(open(info_fn, 'rb'))
        sect_list = [item for item in sect_info.keys()]
        
    else:
        print('** missing sect_info.p **')
        sys.exit()

# gather fields for section extractions
tt0 = time()
ds = nc.Dataset(fn)
ot = ds['ocean_time'][0]
H = ds['h'][:]
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
    ii0, ii1, jj0, jj1, sdir, NX, Lon, Lat = sect_info[sect_name]
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
    C['ot'] = ot
    
    A[sect_name] = C
pickle.dump(A, open(out_fn, 'wb'))
print('Time to extract all section data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

