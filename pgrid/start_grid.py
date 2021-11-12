"""
Code to initialize the creation of a ROMS grid file.

NOTE: the gridname and associated info is all set in LO_user/gfun_user.py.

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to: [lon_rho,lat_rho] as [lon,lat].

Also [plon,plat] is just like [lon_psi, lat_psi] but extended by one on all
directions so that it is box corners around all rho-grid points.

"""
import sys
import pickle
from lo_tools import Lfun

import gfun
import gfun_utility as gfu
Gr =gfun.gstart()
import gfun_user # path added by gfun.gstart()

if Gr['gdir'].is_dir():
    ans = input('Grid ' + Gr['gridname'] + ' exists.  Overwrite? (y/n)')
    if ans == 'y':
        Lfun.make_dir(Gr['gdir'], clean=True)
    else:
        sys.exit()

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = Gr['gdir'] / fn
print(60*'*')
print(str(out_fn).center(60,'-'))

Lfun.make_dir(Gr['gdir'], clean=True)

lon, lat, z, dch = gfun_user.make_initial_info()

# save the output to NetCDF
gfu.make_nc(out_fn, lon, lat, z, dch)

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] / 'choices.p', 'wb'))



