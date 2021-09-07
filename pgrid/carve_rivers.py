"""
This carves the rivers and generates info for ROMS point sources.
"""

from importlib import reload
import gfun; reload(gfun)
Gr =gfun.gstart()
import pfun

import zfun

import os
import shutil
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seawater as sw
import pickle

# load the default choices
dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

#%% get river info

ri_fn = Gr['ri_dir'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

#%% select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_r')
out_fn = Gr['gdir'] + fn_new

#%% get the grid data

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lonu = ds.variables['lon_u'][:]
latu = ds.variables['lat_u'][:]
lonv = ds.variables['lon_v'][:]
latv = ds.variables['lat_v'][:]

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

#%% try carving river channels

m = mask_rho==0

plon_vec = plon[0,:].flatten()
plat_vec = plat[:,0].flatten()

ilon_dict = dict()
ilat_dict = dict()
dir_dict = dict()

good_riv = []

for rn in df.index:
    depth = df.ix[rn, 'depth']
    fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
    
    try:
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        # only include the river if it ENDS in the domain
        if ( x[0] > ax_lims[0] and x[0] < ax_lims[1] and
            y[0] > ax_lims[2] and y[0] < ax_lims[3]):
            good_riv.append(rn)
            print('including ' + rn.title())
            # This unmasks it in the places where the
            # river crosses a tile
            #
            # and we jiggle a little to make sure paths are not blocked
            #nudge_list = [(.003,0), (0,.003), (0,0)]
            nudge_list = [(.0002,0), (0,.0002), (0,0)]
            #nudge_list = [(0,0)]
            x_orig = x.copy()
            y_orig = y.copy()
            for nudge in nudge_list:
                x = x_orig + nudge[0]
                y = y_orig + nudge[1]
                for I in range(len(x)-1):
                    xx = np.linspace(x[I], x[I+1], 100)
                    yy = np.linspace(y[I], y[I+1], 100)
                    ii0, ii1, ifr = zfun.get_interpolant(xx, plon_vec, extrap_nan=True)
                    jj0, jj1, jfr = zfun.get_interpolant(yy, plat_vec, extrap_nan=True)
                    # drop extrapolated points
                    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
                    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
                    # this unmasks points crossed by the river
                    m[jj0, ii0] = False
                    # and this sets the depth in those places, if needed,
                    # "carving the river channel"
                    for ff in range(len(ii0)):
                        JJ = jj0[ff]
                        II = ii0[ff]
                        z[JJ, II] = np.min((z[JJ, II], -depth))
            # this creates information about the index and direction of the river
            ilon_dict[rn] = II
            ilat_dict[rn] = JJ
            # phaseangle is degrees -180:180 with 0 = East
            dist, phaseangle = sw.dist([y[I], y[I+1]], [x[I], x[I+1]])
            dir_dict[rn] = phaseangle
        else:
            print(' >> excluding ' + rn.title())
    except FileNotFoundError:
        pass

#%% figure out river source locations
idir_dict = dict()
isign_dict = dict()
uv_dict = dict()
row_dict_py = dict()
col_dict_py = dict()

ji_dict = dict()
for rn in good_riv:
    ii = ilon_dict[rn]
    jj = ilat_dict[rn]
    ji = np.array([jj,ii])
    ph = dir_dict[rn]
    JI_W = np.array([0,1])
    JI_S = np.array([1,0])
    JI_N = np.array([-1,0])
    JI_E = np.array([0,-1])
    if -45 <= ph and ph <= 45: # River flowing W
        JI = JI_W
    elif 135 <= ph and ph < 45: # S
        JI = JI_S
    elif -135 <= ph and ph < -45: # N
        JI = JI_N
    elif np.abs(ph) > 135: # E
        JI = JI_E
    is_in_grid = True
    while is_in_grid == True:
        # this while loop has the channel keep going until it
        # hits a mask or hits a wall.
        ji = ji + JI
        if (ji[0] == -1) or (ji[1] == -1):
            # we hit the edge of the grid and so move back a step
            ji = ji - JI
            m[ji[0],ji[1]] = True
            ji_dict[rn] = ji
            print(rn + ' (hit a -1 boundary) ' + str(ji))
            is_in_grid = False
        else:
            try:
                # this ensures a minumum depth in the channel
                z[ji[0], ji[1]] = np.min((z[ji[0], ji[1]], -depth))
            except IndexError:
                pass            
            try:
                # we stepped into the masked region
                if m[ji[0],ji[1]] == True:
                    ji_dict[rn] = ji
                    print(rn)
                    is_in_grid = False
            except IndexError:
                # assume we hit the edge of the grid and so move back a step
                ji = ji - JI
                m[ji[0],ji[1]] = True
                ji_dict[rn] = ji
                print(rn + ' (hit boundary) ' + str(ji))
                is_in_grid = False
    if (JI == JI_E).all():
        idir_dict[rn] = 0 # 0 = E/W, 1 = N/S
        isign_dict[rn] = 1
        uv_dict[rn] = 'u'
        row_dict_py[rn] = ji[0]
        col_dict_py[rn] = ji[1]
    elif (JI == JI_W).all():
        idir_dict[rn] = 0
        isign_dict[rn] = -1
        uv_dict[rn] = 'u'
        row_dict_py[rn] = ji[0]
        col_dict_py[rn] = ji[1] - 1
    elif (JI == JI_N).all():
        idir_dict[rn] = 1
        isign_dict[rn] = 1
        uv_dict[rn] = 'v'
        row_dict_py[rn] = ji[0]
        col_dict_py[rn] = ji[1]
    elif (JI == JI_S).all():
        idir_dict[rn] = 1
        isign_dict[rn] = -1
        uv_dict[rn] = 'v'
        row_dict_py[rn] = ji[0] - 1
        col_dict_py[rn] = ji[1]

#%% save all the new info fields to the DataFrame
df = df.ix[good_riv]

df['idir'] = pd.Series(idir_dict)
df['isign'] = pd.Series(isign_dict)
df['uv'] = pd.Series(uv_dict)
df['row_py'] = pd.Series(row_dict_py)
df['col_py'] = pd.Series(col_dict_py)

#%% create the new mask
mask_rho[m == True] = 0
mask_rho[m == False] = 1

#%% plotting

zm = np.ma.masked_where(m, z)

plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
cmap1 = plt.get_cmap(name='rainbow')
cs = ax.pcolormesh(plon, plat,zm,
                   vmin=-200, vmax=20, cmap = cmap1)
fig.colorbar(cs, ax=ax, extend='both')
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(ax_lims)
ax.set_title(Gr['gridname'])

ds.close()

for rn in good_riv:
    fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
    df_tr = pd.read_csv(fn_tr, index_col='ind')
    x = df_tr['lon'].values
    y = df_tr['lat'].values
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[-1], y[-1], '*r')
    
    if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
        ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                latu[row_dict_py[rn], col_dict_py[rn]], '>r')
    elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
        ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                latu[row_dict_py[rn], col_dict_py[rn]], '<r')
    elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
        ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                latv[row_dict_py[rn], col_dict_py[rn]], '^b')
    elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
        ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                latv[row_dict_py[rn], col_dict_py[rn]], 'vb')

plt.show()

#%% Save the output files

print('\nCreating ' + out_fn)
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
shutil.copyfile(in_fn, out_fn)
ds = nc.Dataset(out_fn, 'a')
ds['mask_rho'][:] = mask_rho
ds['h'][:] = -z
ds.close()

out_rfn = Gr['gdir'] + 'river_info.csv'
print('\nCreating ' + out_rfn)
df.to_csv(out_rfn)

# here is how you would read the DataFrame back in
# df1 = pd.read_csv(out_rfn, index_col='rname')

