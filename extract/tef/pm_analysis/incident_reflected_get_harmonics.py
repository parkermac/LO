"""
Get and save all the tidal harmonics for a colleciton of TEF sections.

"""

from pathlib import Path
import sys
import pickle
import xarray as xr
import pandas as pd
import numpy as np
import utide

from lo_tools import Lfun, zfun

gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun

def get_hm(t_days, v, lat, epoch=Lfun.modtime0):
    # function to get harmonics for variable v
    hm = utide.solve(t_days, v, lat=lat, epoch=epoch,
        nodal=True, trend=False)
    return hm
    # hm.aux.freq has units cyles/hour
    # so for f = hm.aux.frq[hm.name == 'M2'][0] we get
    # 1/f = 12.420601202671868 (hours per cycle)
    # hm.A is amplitude, hm.g is phase (degrees)

# set input directory
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
dates = '2018.01.01_2018.12.31'
in_dir = in_dir0 / ('extractions_' + dates)

# get section info
sect_df = tef_fun.get_sect_df(gridname)

# loop over all sections
sect_list = list(sect_df.index)

testing = True
if testing:
    sect_list = ['wb2'] # wb2 has an island in the middle - good mask test
    
hm_e_dict = dict()
hm_u_dict = dict()
for sect_name in sect_list:
    print((' %s ' % (sect_name)).center(60,'-'))
    in_fn = in_dir / (sect_name + '.nc')
    
    # get section info
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    lat = (y0 + y1)/2
    
    # load extracted fields
    ds = xr.open_dataset(in_fn, decode_times=False)
    h = ds['h'].values
    q = ds['q'].values
    DA = ds['DA'].values
    DA0 = ds['DA0'].values
    ot = ds['ocean_time'].values # time in seconds since Lfun.modtime0 (1/1/1970)
    td = ot/86400 # time in days
    zeta = ds['zeta'][:]
    # mask: True for water points on section
    mask = ~np.isnan(q[0,0,:])
    # surface height
    eta = zeta[:,mask].mean(axis=1).data
    # remove low-passed SSH
    eta = eta - zfun.lowpass(eta, f='godin', nanpad=True)
        
    # load extracted fields
    ds = nc.Dataset(in_fn)
    h = ds['h'].values
    q = ds['q'].values
    DA = ds['DA'].values
    DA0 = ds['DA0'][:]
    ot = ds['ocean_time'][:].data
    ot_days = ot/86400
    zeta = ds['zeta'][:]
    NT, NZ, NX = q.shape

    # mask: True for water points on section
    mask = ~np.isnan(q[0,0,:]).data # use q to find the correct mask
    H = h[mask].mean()
    A0 = (DA0.data[:,mask]).sum()

    # velocity
    Q = q.data[:,:,mask].sum(axis=2).sum(axis=1)
    A = DA.data[:,:,mask].sum(axis=2).sum(axis=1)
    u = Q/A

    # surface height
    eta = zeta[:,mask].mean(axis=1).data
    
    # remove low-passed signal
    eta = eta - zfun.lowpass(eta, f='godin', nanpad=False)
    u = u - zfun.lowpass(u, f='godin', nanpad=False)
    # this leaves zeros as the padding on the ends
    
    # calculate harmonics
    hm_e = get_hm(ot_days, eta, lat)
    hm_u = get_hm(ot_days, u, lat)
    
    # and add H, A0, and F for later use
    hm_e['H'] = H
    hm_u['H'] = H
    hm_e['A0'] = A0
    hm_u['A0'] = A0
    rho = 1025
    g = 9.8
    F = rho * g * (Q * eta).mean()
    hm_e['F'] = F
    hm_u['F'] = F
    
    # gather results
    hm_e_dict[sect_name] = hm_e
    hm_u_dict[sect_name] = hm_u
        
if not testing:
    out_dir = in_dir0 / ('harmonics_' + dates)
    Lfun.make_dir(out_dir, clean=True)
    pickle.dump(hm_e_dict, open(out_dir / 'hm_e_dict.p', 'wb'))
    pickle.dump(hm_u_dict, open(out_dir / 'hm_u_dict.p', 'wb'))

