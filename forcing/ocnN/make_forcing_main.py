"""
This makes the ocn forcing files for a nested run, interpolating to get
the fields from history files for another run.

Designed to run only as backfill

Testing:

run make_forcing_main.py -g hc0 -t v0 -r backfill -s continuation -d 2019.07.04 -f ocnN -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
from scipy.spatial import cKDTree
import numpy as np
import matplotlib.pyplot as plt

from lo_tools import Lfun, zfun, zrfun, Ofun_nc
from lo_tools import plotting_functions as pfun

# defaults
verbose = False
if Ldir['testing']:
    verbose = True

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# where to find the files to interpolate from
in_dir = Ldir['parent'] / 'LiveOcean_roms' / 'output' / 'cas6_v3_lo8b' / ('f' + Ldir['date_string'])

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

h_list = sorted(in_dir.glob('ocean_his_*'))

G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
#G0 = zrfun.get_basic_info(h_list[0], only_G=True)

# use nearest neighbor interpolation to fill a field
ds = xr.open_dataset(h_list[0])
x = ds.lon_rho.values
y = ds.lat_rho.values
m = ds.mask_rho.values # 1=water
s = ds.salt[0,-1,:,:].values 
ds.close()

xygood = np.array((x[m==1],y[m==1])).T
xybad = np.array((x[m==0],y[m==0])).T
xyT_rho = cKDTree(xygood)

ss = s.copy()
ss[m==0] = s[m==1][xyT_rho.query(xybad)[1]]

if Ldir['testing']:
    plon, plat = pfun.get_plon_plat(x, y)
    
    plt.close('all')
    pfun.start_plot(figsize=(12,12))
    fig = plt.figure()
    
    ax = fig.add_subplot(121)
    ax.pcolormesh(plon,plat,s)
    pfun.dar(ax)
    
    ax = fig.add_subplot(122)
    ax.pcolormesh(plon,plat,ss)
    pfun.dar(ax)
    
    plt.show()
    pfun.end_plot()

# # and interpolate to ROMS format
# # get grid and S info
# S_fn = Ldir['grid'] / 'S_COORDINATE_INFO.csv'
# S_info_dict = pd.read_csv(S_fn, index_col='ITEMS').to_dict()['VALUES']
# S = zrfun.get_S(S_info_dict)
# # Write to ROMS forcing files
# Ofun_nc.make_clm_file(Ldir, out_dir, h_out_dir, c_dict, dt_list, S, G)
    
    
# Ofun_nc.make_ini_file(out_dir)
# Ofun_nc.make_bry_file(out_dir)

if False:
    # check results
    nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
    result_dict['result'] = 'success'
    for fn in nc_list:
        if (out_dir / fn).is_file():
            pass
        else:
           result_dict['result'] = 'fail'

    # *******************************************************

    result_dict['end_dt'] = datetime.now()
    ffun.finale(Ldir, result_dict)
