"""
Tool for editing the mask, and the depth, of a grid file.  During
depth editing it sets the depth to a constant "dval" set early
in the code, which can also be set as a kwarg.

By using imshow() this is MUCH faster than anything I achieved
using pcolormesh().  E.g with an 800x500 grid it was still
pleasant to use, whereas the old version was unworkable.

"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dval', default=5, type=float)
args = parser.parse_args()

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import sys
import pandas as pd

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import gfun
import gfun_utility as gfu

Gr =gfun.gstart()

Ldir = Lfun.Lstart()

# set the depth to impose during Depth Editing
dval = args.dval # m (positive down)

# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_m')

# get fields
ds = xr.open_dataset(in_fn)
H = ds.h.values
lon = ds.lon_rho.values
lat = ds.lat_rho.values
mask_rho = ds.mask_rho.values
ds.close()

"""
The code works by changing values in an array "hh" which stores
the modified depth with nan's for the modified mask.  When needed
it uses the unmodified depth "h" (which has no nan's) to restore
back to the original depth.

At the end it uses hh to update the depth and mask for the next
NetCDF save of the file.
"""

# flip to work with imshow
h = np.flipud(H)
m = np.flipud(mask_rho)
# mask_rho:
# 1 = water
# 0 = land
hh = h.copy()
hh[m==0] = np.nan
NR, NC = hh.shape

# ----------------- padding ------------------------------------
# Add padding around edges to facilitate editing out areas that
# extend to the boundary.  This is removed at the end before writing to NetCDF.
hpad = 12
NRp = NR + 2*hpad
NCp = NC + 2*hpad

hp = np.zeros((NRp,NCp))
mp = np.zeros((NRp,NCp))
hhp = np.nan * np.ones((NRp,NCp))

hp[hpad:-hpad, hpad:-hpad] = h
hhp[hpad:-hpad, hpad:-hpad] = hh
mp[hpad:-hpad, hpad:-hpad] = m

h = hp.copy()
hh = hhp.copy()
m = mp.copy()

lon = lon[0,:]
lat = lat[:,0]

dlon = lon[1] - lon[0]
dlat = lat[1] - lat[0]

lonp = np.nan * np.ones(NCp)
latp = np.nan * np.ones(NRp)
lonp[hpad:-hpad] = lon
latp[hpad:-hpad] = lat
lonp[:hpad] = np.arange(lon[0] - hpad*dlon, lon[0], dlon)
lonp[-hpad:] = np.arange(lon[-1]+dlon, lon[-1] + (hpad+1)*dlon, dlon)
latp[:hpad] = np.arange(lat[0] - hpad*dlat, lat[0], dlat)
latp[-hpad:] = np.arange(lat[-1]+dlat, lat[-1] + (hpad+1)*dlat, dlat)

lon = lonp.copy()
lat = latp.copy()
# ----------------- padding ------------------------------------

# PLOTTING

# set up the axes
plt.close('all')
fig = plt.figure(figsize=(20,13)) # (13,8) is good for my laptop
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # map
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1) # buttons

# initialize the data plot
# try a segmented colormap
from matplotlib import colors
# make a color map of fixed colors
cmap = colors.ListedColormap(['coral', 'orange', 'lightgreen', 'seagreen',
                              'lightblue', 'cornflowerblue', 'mediumslateblue', 'blueviolet'])
# bounds=[-10, -5, 0, 5, 10, 20, 100, 200, 4000]
bounds=[-5, 0, 1, 2, 3, 5, 10, 20, 4000]
norm = colors.BoundaryNorm(bounds, cmap.N)
# tell imshow about color map so that only set colors are used
cs = ax1.imshow(h, interpolation='nearest', cmap=cmap, norm=norm)
fig.colorbar(cs, ax=ax1)
aa = ax1.axis()

# add the coastline
clon, clat = pfun.get_coast()
# cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:], show_warnings=False)
# cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0], show_warnings=False)
cx0, cx1, cxf = zfun.get_interpolant(clon, lon, show_warnings=False)
cy0, cy1, cyf = zfun.get_interpolant(clat, lat, show_warnings=False)
ax1.plot(cx0 + cxf, NRp - (cy0 + cyf) - 1, '-k')

def edit_mask_river_tracks(Gr, NRp, ax, hpad):
    # add river tracks and endpoints for edit_mask.py
    rri_fn = Gr['gdir'] / 'roms_river_info.csv'
    if rri_fn.is_file():
        df = pd.read_csv(rri_fn, index_col='rname')
    else:
        print('Note from edit_mask_river_tracks(): missing roms_river_info.csv')
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # Plot river endpoints, indicating source direction.  The indexing
    # nudges seem a little non-intuitive, but I believe they are correct.
    for rn in df.index:
        yy = NRp - row_dict_py[rn] - 1 - hpad
        xx = col_dict_py[rn] + hpad
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            # River Source on W side of rho cell
            ax.plot(xx+.5, yy, '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            # River Source on E side of rho cell
            ax.plot(xx+.5, yy, '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            # River Source on S side of rho cell
            ax.plot(xx, yy-.5, '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            # River Source on N side of rho cell
            ax.plot(xx, yy-.5, 'vb')    
# add river info if rivers have been carved
edit_mask_river_tracks(Gr, NRp, ax1, hpad)

ax1.axis(aa)

# Create control buttons
# Note: list is organized from bottom to top
blist = ['start', 'pause', 'continueM', 'continueZ',
         'polyToLand', 'polyToWater', 'lineToLand', 'lineToWater', 'startPoly',
         'undo','discard','done']
# nicer names
Blist = ['Start', 'Pause', 'Edit Mask', 'Edit Depth (' + str(dval) + ' m)',
         'Polygon to Land', 'Polygon to Water',
         'Line to Land', 'Line to Water', 'Start Polygon/Line',
         'Undo', 'Discard ALL Edits', 'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from topography in the ginput
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
def get_indices_in_polygon(plon_poly, plat_poly, NRp, NCp):
    # get indices of points inside a polygon
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    # grid centers
    x = np.arange(NCp)
    y = np.arange(NRp)
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
            hh_prev = hh.copy()
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
            # Note the pushing Pause saves the most recent edits and they
            # can't be undone by the Undo button.
            hh_prev = hh.copy()
            flag_continue = False
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='continueM') and not flag_start:
            hh_prev = hh.copy()
            flag_continue = True
            flag_e = 'm'
            ax1.set_title('EDITING Mask')
        elif (bdict[nb]=='continueZ') and not flag_start:
            hh_prev = hh.copy()
            flag_continue = True
            flag_e = 'z'
            ax1.set_title('EDITING Depth')
        elif (bdict[nb]=='startPoly') and not flag_start:
            hh_prev = hh.copy()
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Poly Points')
            remove_poly()
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='polyToLand') and not flag_start:
            hh_prev = hh.copy()
            flag_continue = False
            ax1.set_title('PAUSED: Changed Poly to Land')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NRp, NCp)
            hh[ji_rho_in[:,0], ji_rho_in[:,1]] = np.nan
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='polyToWater') and not flag_start:
            hh_prev = hh.copy()
            flag_continue = False
            ax1.set_title('PAUSED: Changed Poly to Water')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NRp, NCp)
            hh[ji_rho_in[:,0], ji_rho_in[:,1]] = h[ji_rho_in[:,0], ji_rho_in[:,1]]
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='lineToWater') and not flag_start:
            # This unmasks or carves depth in the places where the
            # line crosses a grid cell - ensures continuous path.
            hh_prev = hh.copy()
            flag_continue = False
            x = np.array(plon_poly)
            y = np.array(plat_poly)
            for ii in range(len(x)-1):
                ix0 = int(x[ii])
                ix1 = int(x[ii+1])
                iy0 = int(y[ii])
                iy1 = int(y[ii+1])
                if ix0==ix1 and iy0==iy1:
                    II = np.array(ix0)
                    JJ = np.array(iy0)
                    hh[JJ, II] = np.max(([JJ, II], dval))
                else:
                    II, JJ = zfun.get_stairstep(ix0, ix1, iy0, iy1)
                    for pp in range(len(II)):
                        hh[JJ[pp], II[pp]] = np.max((h[JJ[pp], II[pp]], dval))
                # this sets the depth of points crossed by the line to dval
                # or the original bathymetry - whichever is deeper.
                ax1.set_title('PAUSED: Changed line to Water of depth >= '
                              + str(dval) + ' m')
            cs.set_data(hh)
            remove_poly()
        elif (bdict[nb]=='lineToLand') and not flag_start:
            # This masks in the places where the
            # line crosses a grid cell - ensures continuous path.
            hh_prev = hh.copy()
            flag_continue = False
            x = np.array(plon_poly)
            y = np.array(plat_poly)
            for ii in range(len(x)-1):
                ix0 = int(x[ii])
                ix1 = int(x[ii+1])
                iy0 = int(y[ii])
                iy1 = int(y[ii+1])
                if ix0==ix1 and iy0==iy1:
                    II = np.array(ix0)
                    JJ = np.array(iy0)
                    hh[JJ, II] = np.nan
                else:
                    II, JJ = zfun.get_stairstep(ix0, ix1, iy0, iy1)
                    for pp in range(len(II)):
                        hh[JJ[pp], II[pp]] = np.nan
                ax1.set_title('PAUSED: Changed line to land')
            cs.set_data(hh)
            remove_poly()
            
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE')
            
        elif (bdict[nb]=='undo') and not flag_start:
            flag_continue = False
            hh = hh_prev.copy()
            cs.set_data(hh)
            ax1.set_title('Paused: Undid last edits')
            
        elif (bdict[nb]=='discard') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE - No changes saved')
            print('No change to mask or bathy')
            sys.exit()
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
# (I should make a discard changes button)

# update fields for output

# remove hpad
h = h[hpad:-hpad, hpad:-hpad]
hh = hh[hpad:-hpad, hpad:-hpad]

h = np.flipud(h)
hh = np.flipud(hh)
newmask = np.ones((NR, NC), dtype=float)
newmask[np.isnan(hh)] = 0.

# We do this to maintain the h field (which will be created from hh when
# we write the output below) as a complete array with no nans.
hh[np.isnan(hh)] = h[np.isnan(hh)]

# save new data and mask, if there are any changes
if np.any(mask_rho != newmask) or np.any(h != hh):
    print('Creating ' + str(out_fn))
    # save the updated mask and z
    ds = xr.open_dataset(in_fn)
    ds.update({'mask_rho': (('eta_rho', 'xi_rho'), newmask)})
    ds.update({'h': (('eta_rho', 'xi_rho'), hh)})
    ds.to_netcdf(out_fn)
    ds.close()
else:
    print('No change to mask or bathy')
