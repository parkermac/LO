"""
Tool for interactively exploring the bathymetry of a
grid file.

Can be run in ipython with a user-specified grid file

run grid_explore.py -g ai0

NOTE: in order to look at an older grid file like cas6, just copy
its grid.nc into LO_output/pgrid/[gridname], and then specify it with the
-g flag when you run this program.

New: This will ask if you want to append sections to an existing collection,
or start a new one, or not save anything.

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import sys
import pandas as pd
import pickle

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import gfun
import gfun_utility as gfu

Gr =gfun.gstart()
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='cas6', type=str)
parser.add_argument('-c', '--collection', default='pm_cas6', type=str)
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
# Set the "small" argument to True to work on a laptop
args = parser.parse_args()
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
# select grid file
in_fn = gfun.select_file(Gr)

# organize the collection to write to
if len(args.collection) > 0:
    use_collection = True
    print('Saving sections to collection: ' + args.collection)
    c_dir = Ldir['LOo'] / 'tef_section_collections'
    Lfun.make_dir(c_dir)
    c_fn = c_dir / (args.collection + '.p')
    if c_fn.is_file():
        print(' - adding to existing collection')
        sect_df = pd.read_pickle(c_fn)
    else:
        print(' - creating new collection')
        sect_df = pd.DataFrame(columns=['pname','lon','lat'])
else:
    use_collection = False
    print('Not saving sections')

# get fields
ds = xr.open_dataset(in_fn)
H = ds.h.values
lon = ds.lon_rho.values
lat = ds.lat_rho.values
mask_rho = ds.mask_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
DA = (1/ds.pm.values) * (1/ds.pn.values)
DA[mask_rho==0] = np.nan
H[mask_rho==0] = np.nan
Hm = np.ma.masked_where(mask_rho==0, H)
aa0 = pfun.get_aa(ds)
ds.close()

# get distances
XM, YM = zfun.ll2xy(lon, lat, np.mean(lon[0,:]), np.mean(lat[:,0]))

# flip to work with imshow
h = np.flipud(H)
da = np.flipud(DA)
xm = np.flipud(XM)
ym = np.flipud(YM)
m = np.flipud(mask_rho) # mask_rho: 1 = water, 0 = land
lonvec = lon[0,:] # no need to flip
latvec = np.flipud(lat[:,0])
NR, NC = h.shape

# PLOTTING

# set up the axes
plt.close('all')
if args.small == True:
    fig = plt.figure(figsize=(13,8))
else:
    fig = plt.figure(figsize=(22,12))
    
extra_map = False # set to True to see a second map on the right with
# correct aspect ratio

if extra_map == True:
    ax1 = plt.subplot2grid((1,5), (0,0), colspan=2) # map
    ax2 = plt.subplot2grid((1,5), (0,2), colspan=1) # buttons
    ax3 = plt.subplot2grid((1,5), (0,3), colspan=2) # second map
else:
    ax1 = plt.subplot2grid((1,4), (0,0), colspan=3) # map
    ax2 = plt.subplot2grid((1,4), (0,3), colspan=1) # buttons
    

# initialize the data plot
cmap1 = plt.get_cmap(name='terrain')
tvmin = 0
tvmax = 200
cs = ax1.imshow(h, interpolation='nearest', vmin=tvmin, vmax=tvmax, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
aa = ax1.axis()

# add the coastline
clon, clat = pfun.get_coast()
cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:], show_warnings=False)
cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0], show_warnings=False)
ax1.plot(cx0 + cxf, NR - (cy0 + cyf) - 1, '-k')

# create control buttons
# list is organized from bottom to top
blist = ['start', 'pause','polyInfo', 'lineInfo',
            'tefSave', 'lineSave',
            'startLine', 'done']
# nicer names
Blist = ['Start', 'Pause','Polygon Info', 'Line Info',
            'Save TEF\nSection', 'Save Line',
            '* Start Line *', 'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

if extra_map == True:
    # initialize the lat,lon map plot
    ax3.pcolormesh(plon, plat, Hm,
        vmin=tvmin, vmax=tvmax, cmap = cmap1)
    pfun.add_coast(ax3)
    ax3.axis(aa0)
    pfun.dar(ax3)

def addButtonLabel(ax, plon, plat, nb, lab, tcol='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((plon[0]+pad,plat[nb]+pad),
                      np.diff(plon)[-1]-2*pad, np.diff(plat)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(plon.mean(),nb, lab, fontsize=12,
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
flag_e = 'p' # 'p' for polygon routines
pline = []
plon_poly = []
plat_poly = []
tef_list = []

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
            cs.set_data(h)
            ax1.set_title('PAUSED')
            # reset button colors
            for bnum in bdict.keys():
                if bdict[bnum] == 'start':
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)
                else:
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
        elif (bdict[nb]=='pause') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='startLine') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Poly Points')
            #remove_poly()
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='polyInfo') and not flag_start:
            flag_continue = False
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NR, NC)
            dap = da[ji_rho_in[:,0], ji_rho_in[:,1]]
            dvp = h[ji_rho_in[:,0], ji_rho_in[:,1]] * da[ji_rho_in[:,0], ji_rho_in[:,1]]
            ap = np.nansum(dap)
            vp = np.nansum(dvp)
            hp = vp/ap
            print('- all values are for unmasked area -')
            print('Volume inside polygon = %0.1f km3' % (vp/1e9) )
            print('Area inside polygon = %0.1f km2' % (ap/1e6) )
            print('Mean Depth inside polygon = %0.1f m' % (hp) )
            inp = input('Push Return to continue\n')
            remove_poly()
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='lineSave') and not flag_start:
            flag_continue = False
            pname = input('Name for saved polygon or line: ')
            poutdir = Ldir['data'] / 'section_lines'
            Lfun.make_dir(poutdir)
            lon_poly = lonvec[plon_poly]
            lat_poly = latvec[plat_poly]
            pdict = {'lon_poly': lon_poly, 'lat_poly': lat_poly}
            poutfn = poutdir / (pname.replace(' ','') + '.p')
            pickle.dump(pdict, open(poutfn, 'wb'))
            print(' - saved to ' + str(poutfn))
            #remove_poly()
            ax1.set_title('PAUSED')
            if extra_map == True:
                ax3.plot(lon_poly, lat_poly, '-*k', linewidth=1)
                ax3.text(np.array(lon_poly).mean(),np.array(lat_poly).mean(),
                    pname, ha='center',va='center')
        elif (bdict[nb]=='lineInfo') and not flag_start:
            flag_continue = False
            x = plon_poly
            y = plat_poly
            dist = 0
            for ii in range(len(x)-1):
                dist += np.sqrt( (xm[y[ii+1],x[ii+1]]-xm[y[ii],x[ii]])**2
                    + (ym[y[ii+1],x[ii+1]]-ym[y[ii],x[ii]])**2 )
            print('Line length = %0.1f km' % (dist/1e3))
            inp = input('Push Return to continue\n')
            remove_poly()
            ax1.set_title('PAUSED')
        
        # ************************************************
        elif (bdict[nb]=='tefSave') and not flag_start:
            flag_continue = False
            ax1.set_title('Name TEF section in terminal', fontweight='bold', color='m')
            plt.draw()
            
            pname = input('Name for saved section: ')
            
            lon_poly = lonvec[plon_poly]
            lat_poly = latvec[plat_poly]
            
            if extra_map == True:
                ax3.plot(lon_poly, lat_poly, '-*k', linewidth=1)
                ax3.text(np.array(lon_poly).mean(),np.array(lat_poly).mean(),pname)
            
            for jj in range(len(lon_poly)):
                sect_df.loc[len(sect_df.index),:] = [pname, lon_poly[jj], lat_poly[jj]]
                
            ax1.set_title('PAUSED')
        # ************************************************
        
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE', fontweight='bold', color='b')
            
            sect_df.to_pickle(c_fn)
                    
        else:
            pass
        plt.draw()

    # this code deals with map input, and only responds when
    # we are clicking on the map
    elif flag_continue and not flag_start:
        # we are in the data field
        if flag_e == 'p':
            # this draws a polygon as you click
            plon_poly.append(ix)
            plat_poly.append(iy)
            remove_poly()
            aa = ax1.axis()
            pline = ax1.plot(plon_poly, plat_poly,'*-r')
            ax1.axis(aa)
            ax1.set_title('Adding to Polygon')
        plt.draw()


