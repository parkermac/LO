"""
Module of functions for plotting grids in pgrid.
"""
import numpy as np
import pandas as pd

    
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
    # add river tracks
    gri_fn = Gr['ri_dir'] / 'river_info.csv'
    df = pd.read_csv(gri_fn, index_col='rname')
    # uv_dict = df['uv']
    # row_dict_py = df['row_py']
    # col_dict_py = df['col_py']
    # isign_dict = df['isign']
    # idir_dict = df['idir']
    # # get grids
    # lon_dict, lat_dict, mask_dict = get_grids(ds)
    # lonu = lon_dict['u']
    # latu = lat_dict['u']
    # lonv = lon_dict['v']
    # latv = lat_dict['v']
    # plot river tracks
    for rn in df.index:
        fn_tr = Gr['ri_dir'] / 'tracks' / (rn + '.p')
        try:
            track_df = pd.read_pickle(fn_tr)
        except FileNotFoundError:
            return
        x = track_df['lon'].to_numpy()
        y = track_df['lat'].to_numpy()
        ax.plot(x, y, '-r', linewidth=2)
        ax.plot(x[-1], y[-1], '*r')    
        # if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
        #     ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
        #             latu[row_dict_py[rn], col_dict_py[rn]], '>r')
        # elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
        #     ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
        #             latu[row_dict_py[rn], col_dict_py[rn]], '<r')
        # elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
        #     ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
        #             latv[row_dict_py[rn], col_dict_py[rn]], '^b')
        # elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
        #     ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
        #             latv[row_dict_py[rn], col_dict_py[rn]], 'vb')
                    
def edit_mask_river_tracks(Gr, NR, ax):
    # add river tracks and endpoints for edit_mask.py
    in_rfn = Gr['gdir'] + 'river_info.csv'
    try:
        df = pd.read_csv(in_rfn, index_col='rname')
    except FileNotFoundError:
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # plot river tracks
    for rn in df.index:
        yy = NR - row_dict_py[rn] - 1
        xx = col_dict_py[rn]
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            ax.plot(xx+.5, yy, '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            ax.plot(xx+.5, yy, '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            ax.plot(xx, yy-.5, '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            ax.plot(xx, yy-.5, 'vb')    

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