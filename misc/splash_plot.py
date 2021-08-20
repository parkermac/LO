"""
A pretty plot for talks.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from PyCO2SYS import CO2SYS
import seawater as sw
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from lo_tools import Lfun, zrfun, zfun
from lo_tools import plotting_functions as pfun

from warnings import filterwarnings
filterwarnings('ignore') # skip a warning from PyCO2SYS

Ldir = Lfun.Lstart()

do_topo = True

# model output
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2021.08.21' / 'ocean_his_0021.nc')
T = zrfun.get_basic_info(fn, only_T=True)
ds = xr.open_dataset(fn)
x = ds.lon_psi.values
y = ds.lat_psi.values
th = ds['temp'][0,-1,1:-1,1:-1].values

if do_topo:
    # topography
    tfn = (Ldir['parent'] / 'ptools_data' / 'topo' /
        'srtm15' / 'topo15.nc')
    tds = xr.open_dataset(tfn)
    step = 3
    tx = tds['lon'][::step].values
    ty = tds['lat'][::step].values
    tz = tds['z'][::step,::step].values
    tz[tz<0] = np.nan

def get_arag(ds, fn, aa, nlev):
    G = zrfun.get_basic_info(fn, only_G=True)
    # find indices that encompass region aa
    i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
    i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
    j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
    j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
    px = G['lon_psi'][j0:j1-1, i0:i1-1]
    py = G['lat_psi'][j0:j1-1, i0:i1-1]
    lat = G['lat_rho'][j0:j1,i0:i1] # used in sw.pres
    # first extract needed fields and save in v_dict
    v_dict = {}
    vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
    for cvn in vn_in_list:
        L = ds[cvn][0,nlev,j0:j1,i0:i1].values
        v_dict[cvn] = L
    # ------------- the CO2SYS steps -------------------------
    # create pressure
    Ld = G['h'][j0:j1,i0:i1]
    Lpres = sw.pres(Ld, lat)
    # get in situ temperature from potential temperature
    Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
    # convert from umol/L to umol/kg using in situ dentity
    Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
    Lalkalinity[Lalkalinity < 100] = np.nan
    LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
    LTIC[LTIC < 100] = np.nan
    CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
        Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
    # PH = CO2dict['pHout']
    # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
    ARAG = CO2dict['OmegaARout']
    ARAG = ARAG.reshape((v_dict['salt'].shape))
    ARAG = ARAG[1:-1, 1:-1]
    return px, py, ARAG

# PLOTTING
plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(12,9))
fig = plt.figure()

# LARGE MAP
ax = fig.add_subplot(121)
cmap = 'RdYlBu_r'
cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
# Inset colorbar
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
fig.colorbar(cs, cax=cbaxes, orientation='vertical')
if do_topo:
    cmap = 'gist_earth'
    cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis([-130, -122, 42, 52])
ax.set_xticks([-129, -127, -125, -123])
ax.set_yticks([42, 44, 46, 48, 50, 52])
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
tstr = T['dt'][0].strftime(Lfun.ds_fmt)
ax.text(.98,.99,'LiveOcean', size=fs*1.2,
     ha='right', va='top', weight='bold', transform=ax.transAxes)
ax.text(.98,.95,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
     ha='right', va='top', weight='bold', transform=ax.transAxes,
     bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

# box for Willapa and Grays Harbor
aa = [-124.6, -123.65, 46, 47.2]
nlev = 0
# draw box on the large map
pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)

# SMALL. MAP
ax = fig.add_subplot(122)
px, py, ARAG = get_arag(ds, fn, aa, nlev)
cs = ax.pcolormesh(px,py,ARAG, cmap='coolwarm_r', vmin=0, vmax=3)
# Inset colorbar
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
fig.colorbar(cs, cax=cbaxes, orientation='vertical')
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(aa)
ax.set_xticks([-124.5, -124])
ax.set_yticks([46, 47])
ax.set_xlabel('Longitude')
ax.text(.98,.99,'Bottom water\nAragonite\nSaturation\nState',
     ha='right', va='top', weight='bold', transform=ax.transAxes)

fig.tight_layout()

if True:
    out_dir = Ldir['LOo'] / 'misc'
    Lfun.make_dir(out_dir)
    plt.savefig(out_dir / 'splash.png', transparent=True)

plt.show()
pfun.end_plot()

ds.close()