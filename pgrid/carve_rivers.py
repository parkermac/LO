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
gri_fn = Gr['ri_dir0'] / dch['ctag'] / 'river_info.csv'
ri_df = pd.read_csv(gri_fn, index_col='rname')

# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_r')

# get the grid data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = ds.mask_rho.values
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:]
Y = lat[:,0]

def in_domain(x, y, X, Y):
    # Utility function to make sure that a point (x, y) is
    # in a domain specified by vectors X and Y.
    # We actually require the point to be 'pad' in from the edge.
    pad = 1
    if x>=X[0+pad] and x<=X[-1-pad] and y>=Y[0+pad] and y<=Y[-1-pad]:
        return True
    else:
        return False

# prepare things for looping over the rivers
ilon_dict = dict()
ilat_dict = dict()
dir_dict = dict()
good_riv = []
roms_info_df = pd.DataFrame()

# loop over rivers
for rn in ri_df.index:
    if rn not in dch['excluded_rivers']:
        depth = ri_df.loc[rn, 'depth']
        try:
            fn_tr = Gr['ri_dir0'] / dch['ctag'] / 'tracks' / (rn + '.p')
            track_df = pd.read_pickle(fn_tr)
            # NOTE: tracks go from ocean to land
            x = track_df['lon'].to_numpy()
            y = track_df['lat'].to_numpy()
        
            if in_domain(x[0], y[0], X, Y):
                # we only consider rivers whose mouth is in the domain
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
                        mask_rho[JJ, II] = 1
                    
                # form unique entries of lists
                JI_ulist = []
                for pp in range(len(II_list)):
                    ji_tup = (JJ_list[pp], II_list[pp])
                    if ji_tup not in JI_ulist:
                        JI_ulist.append(ji_tup)
        
                # Then use the last two point to decide on the river direction.
                # dx and dy are in: [-1,0,1] and are interpreted as, e.g.
                # dx = 1 means that the source is on the E side of the rho cell
                # ( column +1 is the step we took to get to the last river channel rho cell).
                if len(JI_ulist) >1:
                    # simple if we have at least two points in JI_ulist
                    dx = JI_ulist[-1][1] - JI_ulist[-2][1] 
                    dy = JI_ulist[-1][0] - JI_ulist[-2][0]
                else:
                    # more tedious if we have ony one point in JI_ulist
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
        
                # Write info for ROMS based on dx and dy.
                #
                # NOTE: the row_py and col_py that we write to river_info.csv use
                # python indexing conventions, meaning that the rho, u, and v grids
                # start at [0,0].  Hence when we plot things, as in plot_grid.py, we
                # should be able to use [row_py, col_py] directly with the u and v grids.
                #
                # However, when we later create the source positions for running ROMS,
                # e.g. in forcing/riv0/make_forcing_main.py => rivers.nc, we convert these
                # values to ROMS indexing conventions, meaning:
                # - add 1 to row of N/S sources
                # - add 1 to column of E/W sources.
                xoff = 0
                yoff = 0
                if dx==1 and dy==0:
                    # Source on E side of rho-cell
                    idir = 0 # 0 = E/W, 1 = N/S (redundant with 'uv')
                    isign = -1
                    uv = 'u'
                elif dx==-1 and dy==0:
                    # Source on W side of rho-cell
                    idir = 0
                    isign = 1
                    uv = 'u'
                    xoff = -1
                elif dx==0 and dy==1:
                    # Source on N side of rho-cell
                    idir = 1
                    isign = -1
                    uv = 'v'
                elif dx==0 and dy==-1:
                    # Source on S side of rho-cell
                    idir = 1
                    isign = 1
                    uv = 'v'
                    yoff = -1
                else:
                    print(' Inconsistent last points for %s '.center(60,'*') % (rn))
                roms_info_df.loc[rn, 'row_py'] = JI_ulist[-1][0] + yoff
                roms_info_df.loc[rn, 'col_py'] = JI_ulist[-1][1] + xoff
                roms_info_df.loc[rn, 'idir'] = idir
                roms_info_df.loc[rn, 'isign'] = isign
                roms_info_df.loc[rn, 'uv'] = uv
            
            else:
                print(' >> excluding ' + rn.title())
            
        except FileNotFoundError:
            pass

# save the updated mask and z
ds.update({'mask_rho': (('eta_rho', 'xi_rho'), mask_rho)})
ds.update({'h': (('eta_rho', 'xi_rho'), -z)})
ds.to_netcdf(out_fn)
ds.close()

# save the roms river info
out_rfn = Gr['gdir'] / 'roms_river_info.csv'
print('\nCreating ' + str(out_rfn))
roms_info_df.index.name = 'rname'
roms_info_df.to_csv(out_rfn)
# here is how you would read the DataFrame back in
# df1 = pd.read_csv(out_rfn, index_col='rname')

