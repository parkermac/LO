"""
Test whether I can reconstruct a tidal record using constituents from utide.

RESULT: This works well, but requires care to include Greenwich phase and
nodal corrections, both of which I got as tables from Mike Foreman.

General formula from Hal Mofjeld:

h(t) = f*H*cos[ om(t-t0) + vu - G ]

f is the node factor
vu = V+U is the astronomical argument
t0 is the start of the chosen year.

H and G are the amplitude and phase returned by utide,

I got the f and V+U tables from (in Ldir['data'] / 'tides') from
http://www.pac.dfo-mpo.gc.ca/science/oceans/tidal-marees/index-eng.html
but that link appears to be defunct.

"""

from pathlib import Path
import sys
pth = Path(__file__).absolute().parent.parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# imports
import netCDF4 as nc
import pandas as pd
import numpy as np
import utide
import matplotlib.pyplot as plt
import zfun

# NEW
import pytide

import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun

# set input directory
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
year = 2018
dates = str(year) + '.01.01_' + str(year) + '.12.31'
in_dir = in_dir0 / ('extractions_' + dates)

# choose section
sect_name = 'ai1'
print((' %s ' % (sect_name)).center(60,'-'))
in_fn = in_dir / (sect_name + '.nc')

# get section latitude
sect_df = tef_fun.get_sect_df(gridname)
x0, x1, y0, y1 = sect_df.loc[sect_name,:]
lat = (y0 + y1)/2
    
# load extracted fields
ds = nc.Dataset(in_fn)
h = ds['h'][:]
q = ds['q'][:]
DA = ds['DA'][:]
DA0 = ds['DA0'][:]
ot = ds['ocean_time'][:].data # time in seconds since Lfun.modtime0 (1/1/1970)
td = ot/86400 # time in days
T = td - td[0] # for plotting
zeta = ds['zeta'][:]

# mask: True for water points on section
mask = ~np.isnan(q[0,0,:]).data

# surface height
eta = zeta[:,mask].mean(axis=1).data
# remove low-passed SSH
eta = eta - zfun.lowpass(eta, f='godin', nanpad=True)

# calculate harmonics using utide
hm = utide.solve(td, eta, lat=lat, epoch=Lfun.modtime0,
    nodal=True, trend=False)
'''
NOTES on utide:
    
hm.aux.freq has units cyles/hour
so for f = hm.aux.frq[hm.name == 'M2'][0] we get
1/f = 12.420601202671868 (hours per cycle)
    
hm.A is amplitude
    
hm.g is phase (degrees)

I believe that by setting nodal=True the resulting constituents
have the nodal correction accounted for (i.e. removed) so you
must multiply by the nodal corrections from Foreman to get back
to the original signal.
'''

# Load the node factors and astronomical arguments from Foreman
info_dir = Ldir['data'] / 'tide'

f_fn = info_dir / 'Foreman_node_factors.txt'
f_df = pd.read_csv(f_fn, header=0, sep='\t', index_col='Year')
f_ser = f_df.loc[year]

vu_fn = info_dir / 'Foreman_astronomical_arguments.txt'
vu_df = pd.read_csv(vu_fn, header=0, sep='\t', index_col='Year')
vu_ser = vu_df.loc[year,:]

# Make a list that is just the constituents in the tables from Mike Foreman:
# ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
clist = list(vu_ser.index)

# Use the pytide module to get the same information.
# The tutorial is helpful:
# https://pangeo-pytide.readthedocs.io/en/latest/tutorial.html
wt = pytide.WaveTable(clist)
f_alt, vu_alt = wt.compute_nodal_modulations([Lfun.modtime_to_datetime(ot[0])])

# use utide to reconstruct the tide (just for clist)
ER = utide.reconstruct(td, hm, epoch=Lfun.modtime0, constit=clist)
er_utide = ER.h

# reconstruct the tide by hand
er_hand = 0 * ot
tref = ot[0]
for cons in clist:
    f = hm.aux.frq[hm.name == cons][0] # cycles per hour
    om = 2*np.pi*f/3600 # frequency, radians per sec
    E = hm.A[hm.name == cons][0] # amplitude, meters
    Eg = np.pi * hm.g[hm.name == cons][0] /180 # phase, converted to radians
    if (cons in vu_ser.index) and (cons in f_ser.index):
        if True:
            # use Foreman's tables
            vu = np.pi * vu_ser[cons] / 180
            F = f_ser[cons]
        else:
            # get the info from pytide (close to Foremen, but not identical)
            vu = vu_alt[clist.index(cons)]
            F = f_alt[clist.index(cons)]
    else:
        print('Missing info for ' + cons)
        
    # Use the formula from Mofjeld:
    er_hand += F * E * np.cos(om*(ot - tref) + vu - Eg)
    
plt.close('all')
err = er_hand - er_utide
plt.plot(T, er_utide, '-g', T, er_hand, '-m', T, err, '-k', alpha=.8)
plt.show()
print('Error std = %0.3f' % (err.std()))
