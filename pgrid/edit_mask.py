# -*- coding: utf-8 -*-
"""
Tool for editing the mask, and the depth, of a grid file.  During
depth editing it sets the depth to a constant "dval" set early
in the code.

By using imshow() this is MUCH faster than anything I achieved
using pcolormesh().  E.g with an 800x500 grid it was still
pleasant to use, whereas the old version was unworkable.

Now accepts an optional command line argument to set the
carving depth "dval." 2019.04.19

"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dval', nargs='?', default=5, type=float)
args = parser.parse_args()

import gfun
from importlib import reload
reload(gfun)
Gr =gfun.gstart()
# running gfun.gstart() sets the path to include pfun and zfun
import gfun_utility as gfu
reload(gfu)
import gfun_plotting as gfp
reload(gfp)
import pfun
import zfun
reload(zfun)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import os
import shutil
#import pickle

import Lfun
Ldir = Lfun.Lstart()

# set the depth to impose during Depth Editing
dval = args.dval # m (positive down)

flag_testing = False
if not flag_testing:
    # select grid file
    fn = gfun.select_file(Gr)
    in_fn = Gr['gdir'] + fn
elif flag_testing:
    fn = 'grid_m02_r01_s00_x00.nc'
    in_fn = '/Users/PM5/Documents/ptools_output/pgrid/salish1/' + fn

# create new file name
fn_new = gfun.increment_filename(fn, tag='_m')
out_fn = Gr['gdir'] + fn_new

# get fields
ds = nc.Dataset(in_fn)
H = ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]
DA = (1/ds['pm'][:]) * (1/ds['pn'][:])
DA[mask_rho==0] = np.nan
ds.close()

# flip to work with imshow
h = np.flipud(H)
da = np.flipud(DA)
m = np.flipud(mask_rho)
# mask_rho:
# 1 = water
# 0 = land
hh = h.copy()
hh[m==0] = np.nan

NR, NC = hh.shape

# PLOTTING

# set up the axes
plt.close('all')
fig = plt.figure(figsize=(13,8)) # (13,8) is good for my laptop
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # map
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1) # buttons

#%% initialize the data plot

if False:
    cmap1 = plt.get_cmap(name='rainbow_r') # terrain
    tvmin = -20
    tvmax = 200
    cs = ax1.imshow(h, interpolation='nearest', vmin=tvmin, vmax=tvmax, cmap = cmap1)
    fig.colorbar(cs, ax=ax1, extend='both')
else:
    # try a segmented colormap
    from matplotlib import colors
    # make a color map of fixed colors
    cmap = colors.ListedColormap(['red', 'orange', 'lightgreen', 'green',
                                  'cyan', 'blue', 'violet', 'black'])
    bounds=[-10, -5, 0, 5, 10, 20, 100, 200, 4000]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    # tell imshow about color map so that only set colors are used
    cs = ax1.imshow(h, interpolation='nearest',
                        cmap=cmap, norm=norm)
    fig.colorbar(cs, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
aa = ax1.axis()

# add the coastline
clon, clat = pfun.get_coast()
cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:], extrap_nan=True)
cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0], extrap_nan=True)
ax1.plot(cx0 + cxf, NR - (cy0 + cyf) - 1, '-k')

# add rivers
gfp.edit_mask_river_tracks(Gr, NR, ax1)

ax1.axis(aa)

# create control buttons
# list is organized from bottom to top
blist = ['start', 'pause', 'continueM', 'continueZ',
         'polyToLand', 'polyToWater', 'lineToWater', 'startPoly',
         'done']
# nicer names
Blist = ['Start', 'Pause', 'Edit Mask', 'Edit Depth (' + str(dval) + ' m)',
         'Polygon to Land', 'Polygon to Water',
         'Line to Water', 'Start Polygon/Line',
         'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

def addButtonLabel(ax, plon, plat, nb, lab, tcol='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((plon[0]+pad,plat[nb]+pad),
                      np.diff(plon)[-1]-2*pad, np.diff(plat)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(plon.mean(),nb, lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcol)

# label the buttons (numbered bottom to top, 0 to NB-1)
bdict = dict(zip(range(NB), blist))
Bdict = dict(zip(range(NB), Blist))
active_color = 'k'
inactive_color = 'w'
for bnum in bdict.keys():
    if bdict[bnum] == 'start':
        addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
    else:
        addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)

plt.show()

# polygon functions
def get_indices_in_polygon(plon_poly, plat_poly, NR, NC):
    # get indices of points inside a polygon
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    # grid centers
    x = np.arange(NC)
    y = np.arange(NR)
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    M, L = X.shape
    Rlon = X.flatten()
    Rlat = Y.flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    RR = P.contains_points(R) # boolean
    # create arrays of i (column) and j (row) indices
    i_rho = np.arange(L).reshape((1,L)).repeat(M, axis=0)
    j_rho = np.arange(M).reshape((M,1)).repeat(L, axis=1)
    # pack indices that are inside the polygon
    # as a numpy int array, with two columns, packed in order j,i
    ji_rho_in = np.array([j_rho.flatten()[RR], i_rho.flatten()[RR]],
                         dtype=int).T
    return ji_rho_in
def remove_poly():    
    try: # remove old polygon lines if they exist
        pl = pline.pop(0)
        pl.remove()
    except (NameError, IndexError):
        pass

# allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once
flag_e = 'm' # 'm' to edit mask, 'z' to edit z, 'p' for polygon routines
pline = []
plon_poly = []
plat_poly = []

while flag_get_ginput:

    # get ginput, note that you can click with any key
    a = plt.ginput(n=1, timeout=0)
    # returns a list of tuples - of length 1
    b = np.array(a)
    b = np.round(b).astype(int)
    if b.shape != (1,2):
        b = np.array([[-1000, -1000]])
    ix = b[0, 0]
    iy = b[0, 1]

    # this code deals with button input
    if (ix >= offset):
        # were are in the buttons
        nb = iy # button number
        if (bdict[nb]=='start') and flag_start:
            flag_start = False
            flag_continue = False
            cs.set_data(hh)
            ax1.set_title('PAUSED: Initial Mask')
            # reset button colors
            for bnum in bdict.keys():
                if bdict[bnum] == 'start':
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)
                else:
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
        elif (bdict[nb]=='pause') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='continueM') and not flag_start:
            flag_continue = True
            flag_e = 'm'
            ax1.set_title('EDITING Mask')
        elif (bdict[nb]=='continueZ') and not flag_start:
            flag_continue = True
            flag_e = 'z'
            ax1.set_title('EDITING Depth')
        elif (bdict[nb]=='startPoly') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Poly Points')
            remove_poly()
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='polyToLand') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED: Changed Poly to Land')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NR, NC)
            hh[ji_rho_in[:,0], ji_rho_in[:,1]] = np.nan
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='polyToWater') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED: Changed Poly to Water')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NR, NC)
            hh[ji_rho_in[:,0], ji_rho_in[:,1]] = h[ji_rho_in[:,0], ji_rho_in[:,1]]
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='lineToWater') and not flag_start:
            flag_continue = False
            x = np.array(plon_poly)
            y = np.array(plat_poly)
            # This unmasks or carves depth in the places where the
            # line crosses a tile
            #
            # and we jiggle a little to make sure paths are not blocked
            # (note these are indices instead of lat,lon degrees)
            nudge_list = [(1,0), (0,1), (0,0)]
            #nudge_list = [(0,0)]
            x_orig = x.copy()
            y_orig = y.copy()
            for nudge in nudge_list:
                x = x_orig + nudge[0]
                y = y_orig + nudge[1]
                for I in range(len(x)-1):
                    xx = np.linspace(x[I], x[I+1], 100)
                    yy = np.linspace(y[I], y[I+1], 100)
                    ii0, ii1, ifr = zfun.get_interpolant(xx, np.arange(NC), extrap_nan=True)
                    jj0, jj1, jfr = zfun.get_interpolant(yy, np.arange(NR), extrap_nan=True)
                    # drop extrapolated points
                    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
                    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
                    # this sets the depth of points crossed by the line to dval
                    ax1.set_title('PAUSED: Changed line to Water of depth '
                                  + str(dval) + ' m')
                    hh[jj0, ii0] = dval
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE')
        else:
            pass
        plt.draw()

    # this code deals with map input, and only responds when
    # we are clicking on the map
    elif flag_continue and not flag_start:
        # we are in the data field
        if flag_e == 'm':
            # this toggles the colors
            if np.isnan(hh[iy, ix]):
                hh[iy, ix] = h[iy, ix]
            else:
                hh[iy, ix] = np.nan
            cs.set_data(hh)
            ax1.set_title('EDITING: ix=' + str(ix) + ' iy=' + str(iy)
                          + ' h=' + str(int(h[iy, ix])) + ' m')
        elif flag_e == 'z':
            # this carves to a specified depth, and removes the mask
            hh[iy, ix] = dval
            cs.set_data(hh)
            ax1.set_title('EDITING: ix=' + str(ix) + ' iy=' + str(iy)
                          + ' h=' + str(int(h[iy, ix])) + ' m')
        elif flag_e == 'p':
            # this draws a polygon as you click
            plon_poly.append(ix)
            plat_poly.append(iy)
            remove_poly()
            aa = ax1.axis()
            pline = ax1.plot(plon_poly, plat_poly,'*-r')
            ax1.axis(aa)
            ax1.set_title('Adding to Polygon')
        plt.draw()

# Save the output file: executes when you push done button
if True:
    # update fields for output
    h = np.flipud(h)
    hh = np.flipud(hh)
    newmask = np.ones((NR, NC), dtype=float)
    newmask[np.isnan(hh)] = 0.
    hh[np.isnan(hh)] = h[np.isnan(hh)]
    # save new data and mask, if needed
    if np.any(mask_rho != newmask) or np.any(h != hh):
        print('Creating ' + out_fn)
        try:
            os.remove(out_fn)
        except OSError:
            pass # assume error was because the file did not exist
        shutil.copyfile(in_fn, out_fn)
        ds = nc.Dataset(out_fn, 'a')
        ds['h'][:] = hh
        ds['mask_rho'][:] = newmask
        ds.close()
    else:
        print('No change to mask')
