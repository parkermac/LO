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

roms_info_df = pd.DataFrame()

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
            # we only consider rivers that END in the domain
            good_riv.append(rn)
            print('including ' + rn.title())
            II_list = []
            JJ_list = []
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
                        II = np.array(ix0)
                        JJ = np.array(iy0)
                        z[JJ, II] = np.min((z[JJ, II], -depth))
                        II_list.append(ix0)
                        JJ_list.append(iy0)
                    else:
                        II, JJ = zfun.get_stairstep(ix0, ix1, iy0, iy1)
                        for pp in range(len(II)):
                            z[JJ[pp], II[pp]] = np.min((z[JJ[pp], II[pp]], -depth))
                        II_list += list(II)
                        JJ_list += list(JJ)
                    m[JJ, II] = False
                    
                    # form unique entries of lists
                    JI_ulist = []
                    for pp in range(len(II_list)):
                        ji_tup = (JJ_list[pp], II_list[pp])
                        if ji_tup not in JI_ulist:
                            JI_ulist.append(ji_tup)
                    
                    # then use the last two point to decide on the river direction
                    if len(JI_ulist) >1:
                        dx = JI_ulist[-1][1] - JI_ulist[-2][1] 
                        dy = JI_ulist[-1][0] - JI_ulist[-2][0]
                    else:
                        dist, ang = sw.dist([y[0],y[-1]], [x[0], x[-1]])
                        if -45<ang and ang<=45:
                            dx=1; dy = 0
                        elif 45<ang and ang<=135:
                            dx=0; dy = 1
                        elif -135<ang and ang<=-45:
                            dx=0; dy = -1
                        elif ang<=-135 or ang>135:
                            dx=-1; dy = 0
                        else:
                            print(' Inconsistent angle for %s '.center(60,'*') % (rn))
                            
                        if dx==1 and dy==0:
                            idir = 0 # 0 = E/W, 1 = N/S
                            isign = -1
                            uv = 'u'
                        elif dx==-1 and dy==0:
                            idir = 0
                            isign = 1
                            uv = 'u'
                        elif dx==0 and dy==1:
                            idir = 1
                            isign = -1
                            uv = 'v'
                        elif dx==0 and dy==-1:
                            idir = 1
                            isign = 1
                            uv = 'v'
                        else:
                            print(' Inconsistent last points for %s '.center(60,'*') % (rn))
                    roms_info_df.loc[rn, 'row_py'] = JI_ulist[-1][0]
                    roms_info_df.loc[rn, 'col_py'] = JI_ulist[-1][1]
                    roms_info_df.loc[rn, 'idir'] = idir
                    roms_info_df.loc[rn, 'isign'] = isign
                    roms_info_df.loc[rn, 'uv'] = uv
        else:
            print(' >> excluding ' + rn.title())
            
    except FileNotFoundError:
        pass

# save the updated mask
mask_rho_new = mask_rho.copy()
mask_rho_new[m == False] = 1
ds.update({'mask_rho': (('eta_rho', 'xi_rho'), mask_rho_new)})
ds.to_netcdf(out_fn)
ds.close()

# save the river info
out_rfn = Gr['gdir'] / 'river_info.csv'
print('\nCreating ' + str(out_rfn))
roms_info_df.index.name = 'rname'
roms_info_df.to_csv(out_rfn)

# here is how you would read the DataFrame back in
# df1 = pd.read_csv(out_rfn, index_col='rname')

