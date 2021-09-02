"""
More detailed work on getting the incident-reflected calculation
to work.  Focus on section ai1.

"""

from pathlib import Path
import sys
import pickle
import xarray as xr
import pandas as pd
import numpy as np
import utide
import cmath

from lo_tools import Lfun, zfun

gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun

def get_hm(td, E, lat):
    hm = utide.solve(td, E, lat=lat,
        constit=['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2'],
        nodal=False, trend=False, phase='raw')
    return hm

# set input directory
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
dates = '2018.01.01_2018.12.31'
in_dir = in_dir0 / ('extractions_' + dates)

# get section info
sect_df = tef_fun.get_sect_df(gridname)
sect_name = 'ai1'
in_fn = in_dir / (sect_name + '.nc')

# get section info
x0, x1, y0, y1 = sect_df.loc[sect_name,:]
lat = (y0 + y1)/2
    
# load extracted fields
ds = xr.open_dataset(in_fn, decode_times=False)
h = ds['h'].values
q = ds['q'].values
ot = ds['ocean_time'].values
ot_days = ot/86400
zeta = ds['zeta'].values
DA0 = ds['DA0'].values

# mask: True for water points on section
mask = ~np.isnan(q[0,0,:])
# transport
Q = q[:,:,mask].sum(axis=2).sum(axis=1)
# surface height
eta = zeta[:,mask].mean(axis=1)
# remove low-passed signal
eta = eta - zfun.lowpass(eta, f='godin', nanpad=False)
Q = Q - zfun.lowpass(Q, f='godin', nanpad=False)
# this leaves zeros as the padding on the ends
A0 = (DA0[:,mask]).sum()
H = h[mask].mean()

# calculate harmonics
hme = get_hm(ot_days, eta, lat)
hmq = get_hm(ot_days, Q, lat)

# calculate the tidal energy flux
rho = 1025
g = 9.8
F = rho * g * (Q * eta).mean()

# reconstruct the tidal energy flux from the harmonics
FF = dict()
CE = dict()
CQ = dict()
Fr = 0
for cons in hme.name:
    Ae = hme.A[hme.name == cons][0]
    Aq = hmq.A[hmq.name == cons][0]
    Ge = np.pi * hme.g[hme.name == cons][0] / 180
    Gq = np.pi * hmq.g[hmq.name == cons][0] / 180
    Ce = cmath.rect(Ae, Ge)  # complex amplitude of eta
    Cq = cmath.rect(Aq, Gq)  # complex amplitude of Q
    CE[cons] = Ce
    CQ[cons] = Cq
    ff = rho * g * (Ce.real*Cq.real + Ce.imag*Cq.imag)/2
    Fr += ff
    FF[cons] = ff
print('F = %0.1f, FF = %0.1f' % (F/1e6, Fr/1e6)) # same to within 2%


ce = CE['M2']
cq = CQ['M2']
cf = FF['M2']
cfr = rho * g * (ce.real*cq.real + ce.imag*cq.imag)/2
print('cf = %0.1f, cfr = %0.1f' % (cf/1e6, cfr/1e6)) # identical of course

for fric in [0, .1, .5, 1, 2, 5, 10]:
    alpha = np.sqrt(g/H)/np.sqrt(1 + fric*1j)

    cep = (ce + (cq/A0)/alpha)/2
    cem = (ce - (cq/A0)/alpha)/2

    cqp = A0*cep*alpha
    cqm = A0*cem*alpha
    
    cfp = rho * g * (cep.real*cqp.real + cep.imag*cqp.imag)/2 
    cfm = - rho * g * (cem.real*cqm.real + cem.imag*cqm.imag)/2 
    cfpm = rho*g*(cem.real*cqp.real + cem.imag*cqp.imag)/2 - rho*g*(cep.real*cqm.real + cep.imag*cqm.imag)/2
    cfrr = cfp + cfm
    print('fric = %0.1f, cfp = %0.1f, cfm = %0.1f, cfpm = %0.1f, cfrr = %0.1f'
        % (fric, cfp/1e6, cfm/1e6, cfpm/1e6, cfrr/1e6))




    