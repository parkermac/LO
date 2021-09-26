"""
This is to check on the exact location and sign of river sources for different
instances of pgrid.

RESULT:
in the new hc0 grid the river positions are one step into the river!

"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

# grid made with original ptools/pgrid and rivers from LiveOcean/forcing/riv2
old_fng = Ldir['parent'] / 'LiveOcean_data' / 'grids' / 'cas6' / 'grid.nc'
old_fnr = Ldir['parent'] / 'LiveOcean_output' / 'cas6_v3' / 'f2019.07.04' / 'riv2' / 'rivers.nc'

# grid made with LO/pgrid and rivers from LO/forcing/riv0
new_fng = Ldir['data'] / 'grids' / 'hc0' / 'grid.nc'
new_fnr = Ldir['LOo'] / 'forcing' / 'hc0_v0' / 'f2019.07.04' / 'riv0' / 'rivers.nc'

# get old river name list and etc.
dsr0 = xr.open_dataset(old_fnr)
orn = dsr0['river_name'].values
NR = orn.shape[1]
rn0_list = []
for ii in range(NR):
    a = orn[:,ii]
    if isinstance(a, np.ma.MaskedArray):
        a = a.data
    r = []
    for l in a:
        r.append(l.decode())
    rr = ''.join(r)
    rn0_list.append(rr)
rX0 = dsr0.river_Xposition.values.astype(int)
rY0 = dsr0.river_Eposition.values.astype(int)
rd0 = dsr0.river_direction.values.astype(int) # 0=u-grid, 1=v=grid
dsr0.close()

# get new river name list and etc.
dsr1 = xr.open_dataset(new_fnr)
rn1_list = list(dsr1['river_name'].values)
rX1 = dsr1.river_Xposition.values.astype(int)
rY1 = dsr1.river_Eposition.values.astype(int)
rd1 = dsr1.river_direction.values.astype(int)
dsr1.close()

# grids
def get_grids(fn):
    dsg = xr.open_dataset(fn)
    x = dsg.lon_rho.values
    y = dsg.lat_rho.values
    m = dsg.mask_rho.values
    xp, yp = pfun.get_plon_plat(x, y)
    xu = dsg.lon_u.values
    yu = dsg.lat_u.values
    xv = dsg.lon_v.values
    yv = dsg.lat_v.values
    dsg.close()
    return x,y,m,xp,yp,xu,yu,xv,yv
    
x0,y0,m0,xp0,yp0,xu0,yu0,xv0,yv0 = get_grids(old_fng)
x1,y1,m1,xp1,yp1,xu1,yu1,xv1,yv1 = get_grids(new_fng)

# PLOTTING
ms = 2

plt.close('all')
fig = plt.figure(figsize=(20,10))

ax = fig.add_subplot(121)
ax.plot(x0[m0==1],y0[m0==1],'oc', ms=ms)
# the "-1" in the indexing is to convert from ROMS 1-based numbering to python 0-based
ax.plot(xu0[rY0[rd0==0],rX0[rd0==0]-1], yu0[rY0[rd0==0],rX0[rd0==0]-1],'>r', ms=ms)
ax.plot(xv0[rY0[rd0==1]-1,rX0[rd0==1]], yv0[rY0[rd0==1]-1,rX0[rd0==1]],'^b', ms=ms)

ax = fig.add_subplot(122)
ax.plot(x1[m1==1],y1[m1==1],'oc', ms=ms)
ax.plot(xu1[rY1[rd1==0],rX1[rd1==0]-1], yu1[rY1[rd1==0],rX1[rd1==0]-1],'>r', ms=ms)
ax.plot(xv1[rY1[rd1==1]-1,rX1[rd1==1]], yv1[rY1[rd1==1]-1,rX1[rd1==1]],'^b', ms=ms)

plt.show()



