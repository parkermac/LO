"""
Code to initialize the creation of a ROMS grid file.

NOTE: the gridname and associated info is all set in LO_user/gfun_user.py.

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to: [lon_rho,lat_rho] as [lon,lat].

Also [plon,plat] is just like [lon_psi, lat_psi] but extended by one on all
directions so that it is box corners around all rho-grid points.

"""
import numpy as np
import pickle
from lo_tools import Lfun, zfun

import gfun
import gfun_utility as gfu
import gfun_user

testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfu)
    reload(gfun_user)
    
Gr =gfun.gstart()
Lfun.make_dir(Gr['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = Gr['gdir'] / fn
print(60*'*')
print(str(out_fn).center(60,'-'))

dch = gfun.default_choices(Gr)
lon, lat, z, dch = gfun_user.make_initial_info(dch)

# save the output to NetCDF
gfu.make_nc(out_fn, lon, lat, z, dch)

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] / 'choices.p', 'wb'))



