"""
This carves the rivers and generates info for ROMS point sources.
"""
import gfun
Gr =gfun.gstart()

from lo_tools import zfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seawater as sw
import pickle

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# get river info
gri_fn = Gr['ri_dir'] / 'river_info.csv'
df = pd.read_csv(gri_fn, index_col='rname')

#% select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] / fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_r')
out_fn = Gr['gdir'] / fn_new

# get the grid data

ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = ds.mask_rho.values
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:]
Y = lat[:,0]
lonu = ds.lon_u.values
latu = ds.lat_u.values
lonv = ds.lon_v.values
latv = ds.lat_v.values

plon, plat = pfun.get_plon_plat(lon, lat)

m = mask_rho==0 # True where masked

ilon_dict = dict()
ilat_dict = dict()
dir_dict = dict()

good_riv = []

def in_domain(xx, yy, XX, YY):
    if xx>XX[0] and xx<XX[-1] and yy>YY[0] and yy<YY[-1]:
        return True
    else:
        return False
        
for rn in df.index:
    depth = df.loc[rn, 'depth']
    
    try:
        fn_tr = Gr['ri_dir'] / 'tracks' / (rn + '.p')
        track_df = pd.read_pickle(fn_tr)
        x = track_df['lon'].to_numpy()
        y = track_df['lat'].to_numpy()
        
        if in_domain(x[0], y[0], X, Y):
            good_riv.append(rn)
            print('including ' + rn.title())
            
        for ii in range(len(x)-1):
            x0 = x[ii]
            x1 = x[ii+1]
            y0 = y[ii]
            y1 = y[ii+1]

            if in_domain(x0, y0, X, Y) and in_domain(x1, y1, X, Y):
                ix0 = zfun.find_nearest_ind(X,x0)
                ix1 = zfun.find_nearest_ind(X,x1)
                iy0 = zfun.find_nearest_ind(Y,y0)
                iy1 = zfun.find_nearest_ind(Y,y1)
                if ix0==ix1 and iy0==iy1:
                    II = ix0
                    JJ = iy0
                else:
                    II, JJ = zfun.get_stairstep(ix0, ix1, iy0, iy1)
                m[JJ, II] = False
                #z[JJ, II] = np.min((z[JJ, II], -depth))

            # # this creates information about the index and direction of the river
            # ilon_dict[rn] = II
            # ilat_dict[rn] = JJ
            # # phaseangle is degrees -180:180 with 0 = East
            # dist, phaseangle = sw.dist([y[I], y[I+1]], [x[I], x[I+1]])
            # dir_dict[rn] = phaseangle
        else:
            print(' >> excluding ' + rn.title())
            
    except FileNotFoundError:
        pass

mask_rho_new = mask_rho.copy()
mask_rho_new[m == False] = 1
ds.update({'mask_rho': (('eta_rho', 'xi_rho'), mask_rho_new)})
ds.to_netcdf(out_fn)
ds.close()

# #%% figure out river source locations
# idir_dict = dict()
# isign_dict = dict()
# uv_dict = dict()
# row_dict_py = dict()
# col_dict_py = dict()
#
# ji_dict = dict()
# for rn in good_riv:
#     ii = ilon_dict[rn]
#     jj = ilat_dict[rn]
#     ji = np.array([jj,ii])
#     ph = dir_dict[rn]
#     JI_W = np.array([0,1])
#     JI_S = np.array([1,0])
#     JI_N = np.array([-1,0])
#     JI_E = np.array([0,-1])
#     if -45 <= ph and ph <= 45: # River flowing W
#         JI = JI_W
#     elif 135 <= ph and ph < 45: # S
#         JI = JI_S
#     elif -135 <= ph and ph < -45: # N
#         JI = JI_N
#     elif np.abs(ph) > 135: # E
#         JI = JI_E
#     is_in_grid = True
#     while is_in_grid == True:
#         # this while loop has the channel keep going until it
#         # hits a mask or hits a wall.
#         ji = ji + JI
#         if (ji[0] == -1) or (ji[1] == -1):
#             # we hit the edge of the grid and so move back a step
#             ji = ji - JI
#             m[ji[0],ji[1]] = True
#             ji_dict[rn] = ji
#             print(rn + ' (hit a -1 boundary) ' + str(ji))
#             is_in_grid = False
#         else:
#             try:
#                 # this ensures a minumum depth in the channel
#                 z[ji[0], ji[1]] = np.min((z[ji[0], ji[1]], -depth))
#             except IndexError:
#                 pass
#             try:
#                 # we stepped into the masked region
#                 if m[ji[0],ji[1]] == True:
#                     ji_dict[rn] = ji
#                     print(rn)
#                     is_in_grid = False
#             except IndexError:
#                 # assume we hit the edge of the grid and so move back a step
#                 ji = ji - JI
#                 m[ji[0],ji[1]] = True
#                 ji_dict[rn] = ji
#                 print(rn + ' (hit boundary) ' + str(ji))
#                 is_in_grid = False
#     if (JI == JI_E).all():
#         idir_dict[rn] = 0 # 0 = E/W, 1 = N/S
#         isign_dict[rn] = 1
#         uv_dict[rn] = 'u'
#         row_dict_py[rn] = ji[0]
#         col_dict_py[rn] = ji[1]
#     elif (JI == JI_W).all():
#         idir_dict[rn] = 0
#         isign_dict[rn] = -1
#         uv_dict[rn] = 'u'
#         row_dict_py[rn] = ji[0]
#         col_dict_py[rn] = ji[1] - 1
#     elif (JI == JI_N).all():
#         idir_dict[rn] = 1
#         isign_dict[rn] = 1
#         uv_dict[rn] = 'v'
#         row_dict_py[rn] = ji[0]
#         col_dict_py[rn] = ji[1]
#     elif (JI == JI_S).all():
#         idir_dict[rn] = 1
#         isign_dict[rn] = -1
#         uv_dict[rn] = 'v'
#         row_dict_py[rn] = ji[0] - 1
#         col_dict_py[rn] = ji[1]
#
# #%% save all the new info fields to the DataFrame
# df = df.ix[good_riv]
#
# df['idir'] = pd.Series(idir_dict)
# df['isign'] = pd.Series(isign_dict)
# df['uv'] = pd.Series(uv_dict)
# df['row_py'] = pd.Series(row_dict_py)
# df['col_py'] = pd.Series(col_dict_py)
#
# #%% create the new mask
# mask_rho[m == True] = 0
# mask_rho[m == False] = 1
#
# #%% plotting
#
# zm = np.ma.masked_where(m, z)
#
# plt.close()
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111)
# cmap1 = plt.get_cmap(name='rainbow')
# cs = ax.pcolormesh(plon, plat,zm,
#                    vmin=-200, vmax=20, cmap = cmap1)
# fig.colorbar(cs, ax=ax, extend='both')
# pfun.add_coast(ax)
# pfun.dar(ax)
# ax.axis(ax_lims)
# ax.set_title(Gr['gridname'])
#
# ds.close()
#
# for rn in good_riv:
#     fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
#     df_tr = pd.read_csv(fn_tr, index_col='ind')
#     x = df_tr['lon'].values
#     y = df_tr['lat'].values
#     ax.plot(x, y, '-r', linewidth=2)
#     ax.plot(x[-1], y[-1], '*r')
#
#     if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
#         ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
#                 latu[row_dict_py[rn], col_dict_py[rn]], '>r')
#     elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
#         ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
#                 latu[row_dict_py[rn], col_dict_py[rn]], '<r')
#     elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
#         ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
#                 latv[row_dict_py[rn], col_dict_py[rn]], '^b')
#     elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
#         ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
#                 latv[row_dict_py[rn], col_dict_py[rn]], 'vb')
#
# plt.show()
#
# #%% Save the output files
#
# print('\nCreating ' + out_fn)
# try:
#     os.remove(out_fn)
# except OSError:
#     pass # assume error was because the file did not exist
# shutil.copyfile(in_fn, out_fn)
# ds = nc.Dataset(out_fn, 'a')
# ds['mask_rho.values = mask_rho
# ds['h.values = -z
# ds.close()
#
# out_rfn = Gr['gdir'] + 'river_info.csv'
# print('\nCreating ' + out_rfn)
# df.to_csv(out_rfn)
#
# # here is how you would read the DataFrame back in
# # df1 = pd.read_csv(out_rfn, index_col='rname')

