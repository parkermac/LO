"""
Module of functions for plotting grids in pgrid.
"""
import numpy as np
import pandas as pd
import pickle

    
def get_grids(ds):
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
    return (lon_dict, lat_dict, mask_dict)
    
def add_river_tracks(Gr, ds, ax):
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    lon_u = ds.lon_u.values
    lat_u = ds.lat_u.values
    lon_v = ds.lon_v.values
    lat_v = ds.lat_v.values
    
    rri_fn = Gr['gdir'] / 'roms_river_info.csv'
    rri_df = pd.read_csv(rri_fn, index_col='rname')
    
    # load the default choices
    dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

    for rn in rri_df.index:
        # These are indices (python, zero-based) into either the
        # u or v grids.
        ii = int(rri_df.loc[rn,'col_py'])
        jj = int(rri_df.loc[rn,'row_py'])
        
        uv = rri_df.loc[rn,'uv']
        isign = rri_df.loc[rn,'isign']
        idir = rri_df.loc[rn,'idir']
        
        if uv == 'u' and isign == 1:
            # River source on W side of rho cell
            ax.plot(lon_u[jj,ii], lat_u[jj,ii],'>r')
            ax.plot(lon[jj,ii+1], lat[jj,ii+1],'oc')
        if uv == 'u' and isign == -1:
            # River source on E side of rho cell
            ax.plot(lon_u[jj,ii], lat_u[jj,ii],'<r')
            ax.plot(lon[jj,ii], lat[jj,ii],'oc')
        if uv == 'v' and isign == 1:
            # River source on S side of rho cell
            ax.plot(lon_v[jj,ii], lat_v[jj,ii],'^b')
            ax.plot(lon[jj+1,ii], lat[jj+1,ii],'oc')
        if uv == 'v' and isign == -1:
            # River source on N side of rho cell
            ax.plot(lon_v[jj,ii], lat_v[jj,ii],'vb')
            ax.plot(lon[jj,ii], lat[jj,ii],'oc')
            
        # fn_tr = Gr['ri_dir0'] / dch['ctag'] / 'tracks' / (rn + '.p')
        # try:
        #     track_df = pd.read_pickle(fn_tr)
        # except FileNotFoundError:
        #     return
        # x = track_df['lon'].to_numpy()
        # y = track_df['lat'].to_numpy()
        # ax.plot(x, y, '-r', linewidth=1, alpha=.3)
        # ax.plot(x[-1], y[-1], '*r', alpha=.3)
                    
def show_z_info(zm, ax):
    # find the max value of z (DEBUGGING)
    (rowmax, colmax) = np.unravel_index(np.argmax(zm), zm.shape)
    zmax = zm[rowmax, colmax]
    print('Max z = ' + str(zmax))
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    ax.plot(lon_rho[rowmax, colmax], lat_rho[rowmax, colmax], '*m', markersize=20)

def show_grids(ds, ax):
    lon_dict, lat_dict, mask_dict = get_grids(ds)
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))