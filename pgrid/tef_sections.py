"""
Tool for interactively defining TEF sections.

Can be run in ipython with a user-specified grid file

run tef_sections.py -g ai0

NOTE: in order to look at an older grid file like cas6, just copy
its grid.nc into LO_output/pgrid/[gridname], and then specify it with the
-g flag when you run this program.

New: This is now designed to be a more useful tool for generating a collection
of TEF sections. Using the -c argument you can use an existing collection or
start a new one. Buttons allow you save a new polygon to a named TEF section,
and we use a TextBox widget to allow you to name the section. There is also
a button to allow you delete any of the sections that were in the existing
collection. I don't yet have a way to delete the sections created in a session
(but you can delete them in the next session).

The TEF sections are saved in a pickled pandas DataFrame, saving just the indices
into the grid for all the points defining the section (one row per point).
It assumes that these will be used for a new TEF section extractor that allows
for diagonal and multi-segment lines.

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import sys
import pandas as pd
import pickle
from matplotlib.widgets import TextBox, Button

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import gfun
import gfun_utility as gfu

Gr =gfun.gstart()
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
# name the grid and the collection
parser.add_argument('-g', '--gridname', default='ae0', type=str)
parser.add_argument('-c', '--collection', default='', type=str)
# set clobber to True to start a new collection even if it exists
parser.add_argument('-clobber', default=False, type=Lfun.boolean_string)
# set small to True to work on a laptop
parser.add_argument('-small', default=True, type=Lfun.boolean_string)
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
    temp_dir = c_dir / 'temp' # the temp dir is for new sections
    Lfun.make_dir(c_dir)
    Lfun.make_dir(temp_dir, clean=True)
    c_fn = c_dir / (args.collection + '.p')
    if c_fn.is_file() and (args.clobber == False):
        print(' - adding to existing collection')
        sect_df = pd.read_pickle(c_fn)
    else:
        print(' - creating new collection')
        c_fn.unlink(missing_ok=True)
        sect_df = pd.DataFrame(columns=['sn','ilon','ilat'])
else:
    use_collection = False
    print('Not saving sections')

# get fields
ds = xr.open_dataset(in_fn)
H = ds.h.values
lon = ds.lon_rho.values
lat = ds.lat_rho.values
mask_rho = ds.mask_rho.values
H[mask_rho==0] = np.nan
ds.close()

# flip to work with imshow
h = np.flipud(H)
NR, NC = h.shape

# PLOTTING

# set up the axes
plt.close('all')
if args.small == True:
    fig = plt.figure(figsize=(13,7.5))
else:
    fig = plt.figure(figsize=(22,12))
    
ax1 = plt.subplot2grid((1,4), (0,0), colspan=3) # map
ax2 = plt.subplot2grid((8,4), (0,3), rowspan=3, colspan=1) # buttons
# box to type strings into
start_ax = plt.subplot2grid((8,4), (5,3))
tb_name_ax = plt.subplot2grid((8,4), (6,3), rowspan=1, colspan=1)
tb_remove_ax = plt.subplot2grid((8,4), (7,3), rowspan=1, colspan=1)
# the widgets
start_button = Button(start_ax, 'Start', image=None, color='0.85', hovercolor='0.95')
tb_name = TextBox(tb_name_ax, 'Name\nSection', initial='', color='.95',
    hovercolor='1', label_pad=0.03)
tb_remove = TextBox(tb_remove_ax, 'Remove\nSection', initial='', color='.95',
    hovercolor='1', label_pad=0.03)
    
# initialize the data plot
cmap1 = plt.get_cmap(name='terrain_r')
tvmin = -20
tvmax = 200
cs = ax1.imshow(h, interpolation='nearest', vmin=tvmin, vmax=tvmax, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
aa = ax1.axis()

# add the coastline
clon, clat = pfun.get_coast()
cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:], show_warnings=False)
cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0], show_warnings=False)
ax1.plot(cx0 + cxf, NR - (cy0 + cyf) - 1, '-k')

def draw_sections(fn):
    # initialize some dicts to allow us to delete lines on the plot
    line_dict = dict()
    text_dict = dict()
    df = pd.read_pickle(fn)
    for sn in df.sn.unique():
        sdf = df[df.sn == sn]
        ilon_vec = sdf.ilon.to_numpy()
        ilat_vec = sdf.ilat.to_numpy()
        line_dict[sn] = ax1.plot(ilon_vec, ilat_vec, '-*b', linewidth=1)
        text_dict[sn] = ax1.text(ilon_vec[0], ilat_vec[0],sn,color='b')
    return line_dict, text_dict


W = dict()
W['flag_get_ginput'] = True # Make False to exit the ginput loop
W['flag_continue'] = False # need to push START to make this True
W['flag_start'] = True # to ensure we only push the start button once
W['flag_e'] = 'p' # 'p' for polygon routines
# pline = []
# ilon_vec = []
# ilat_vec = []

def on_start(junk):
    W['flag_start'] = False
    W['flag_continue'] = False
    cs.set_data(h)
    ax1.set_title('PAUSED')
start_button.on_clicked(on_start)
    
def name_section(sect_name):
    temp_df = pd.DataFrame(columns=['ilon','ilat'])
    temp_df.loc[:,'ilon'] = ilon_vec
    temp_df.loc[:,'ilat'] = ilat_vec
    if len(sect_name) > 0:
        temp_df.to_pickle(temp_dir / (sect_name + '.p'))
    remove_poly(pline)
    ax1.plot(ilon_vec, ilat_vec, '-*k', linewidth=1)
    ax1.text(ilon_vec[0], ilat_vec[0], sect_name)
    tb_name.set_val('')
    ax1.set_title('PAUSED', color='g')
    plt.draw()
tb_name.on_submit(name_section)

def remove_section(x_sect_name):
    dfr = pd.read_pickle(c_fn)
    xx = dfr.sn==x_sect_name
    dfr = dfr.drop(index=dfr.index[xx], axis=1)
    # df = df.reset_index()
    dfr.to_pickle(c_fn)
    try:
        remove_poly(line_dict[x_sect_name])
        text_dict[x_sect_name].remove()
    except:
        pass
    tb_remove.set_val('')
    ax1.set_title('PAUSED', color='g')
    plt.draw()
tb_remove.on_submit(remove_section)

# add existing sections
if use_collection and c_fn.is_file():
    line_dict, text_dict = draw_sections(c_fn)

# create control buttons
# list is organized from bottom to top
blist = ['start', 'startLine', 'done']
# nicer names
Blist = ['Start', '* Start Line *', 'Done']
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
                      fill=True, facecolor=inactive_color, edgecolor='w'))
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
    
def remove_poly(this_line):
    try: # remove old polygon lines if they exist
        pl = this_line.pop(0)
        pl.remove()
    except (KeyError, NameError, IndexError):
        pass

# allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once
flag_e = 'p' # 'p' for polygon routines
pline = []
ilon_vec = []
ilat_vec = []

while flag_get_ginput:

    # get ginput, note that you can click with any key
    a = plt.ginput(n=1, timeout=0)
    # returns a list of tuples - of length 1
    if len(a) == 0:
        # use return to exit line generation
        flag_continue = False
        ax1.set_title('Name TEF section', fontweight='bold', color='m')
        plt.draw()
        
    else:
        b = a[0]
        ix = np.round(b[0]).astype(int)
        iy = np.round(b[1]).astype(int)
    

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
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum],
                        tcol=inactive_color)
                else:
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum],
                        tcol=active_color)
                        
        elif (bdict[nb]=='startLine') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Points (Return to Stop)')
            pline = []
            ilon_vec = []
            ilat_vec = []
        
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE', fontweight='bold', color='b')
            
            if c_fn.is_file() and use_collection:
                sect_df = pd.read_pickle(c_fn)
            else:
                sect_df = pd.DataFrame(columns=['sn','ilon','ilat'])
            
            for fn in temp_dir.glob('*.p'):
                df = pd.read_pickle(fn)
                sn = fn.name.replace('.p','')
                ilon_vec = df.ilon.to_numpy()
                ilat_vec = df.ilat.to_numpy()
                for jj in range(len(ilon_vec)):
                    d= {'sn':sn,'ilon':ilon_vec[jj],'ilat':ilat_vec[jj]}
                    sect_df = sect_df.append(d,ignore_index=True)
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
            ilon_vec.append(ix)
            ilat_vec.append(iy)
            remove_poly(pline)
            aa = ax1.axis()
            pline = ax1.plot(ilon_vec, ilat_vec,'*-r')
            ax1.axis(aa)
            ax1.set_title('Adding Points (Return to Stop)')
        plt.draw()
        
# SPARE CODE that may be useful in the future

# # get fields
# ds = xr.open_dataset(in_fn)
# H = ds.h.values
# lon = ds.lon_rho.values
# lat = ds.lat_rho.values
# mask_rho = ds.mask_rho.values
# plon, plat = pfun.get_plon_plat(lon,lat)
# DA = (1/ds.pm.values) * (1/ds.pn.values)
# DA[mask_rho==0] = np.nan
# H[mask_rho==0] = np.nan
# Hm = np.ma.masked_where(mask_rho==0, H)
# aa0 = pfun.get_aa(ds)
# ds.close()
#
# # get distances
# XM, YM = zfun.ll2xy(lon, lat, np.mean(lon[0,:]), np.mean(lat[:,0]))
#
# # flip to work with imshow
# h = np.flipud(H)
# da = np.flipud(DA)
# xm = np.flipud(XM)
# ym = np.flipud(YM)
# m = np.flipud(mask_rho) # mask_rho: 1 = water, 0 = land
# lonvec = lon[0,:] # no need to flip
# latvec = np.flipud(lat[:,0])
# NR, NC = h.shape

# polygon functions
# def get_indices_in_polygon(ilon_vec, ilat_vec, NR, NC):
#     # get indices of points inside a polygon
#     V = np.ones((len(ilon_vec),2))
#     V[:,0] = ilon_vec
#     V[:,1] = ilat_vec
#     P = mpath.Path(V)
#     # grid centers
#     x = np.arange(NC)
#     y = np.arange(NR)
#     # matrix versions of grids
#     X, Y = np.meshgrid(x,y)
#     M, L = X.shape
#     Rlon = X.flatten()
#     Rlat = Y.flatten()
#     R = np.ones((len(Rlon),2))
#     R[:,0] = Rlon
#     R[:,1] = Rlat
#     RR = P.contains_points(R) # boolean
#     # create arrays of i (column) and j (row) indices
#     i_rho = np.arange(L).reshape((1,L)).repeat(M, axis=0)
#     j_rho = np.arange(M).reshape((M,1)).repeat(L, axis=1)
#     # pack indices that are inside the polygon
#     # as a numpy int array, with two columns, packed in order j,i
#     ji_rho_in = np.array([j_rho.flatten()[RR], i_rho.flatten()[RR]],
#                          dtype=int).T
#     return ji_rho_in

# elif (bdict[nb]=='polyInfo') and not flag_start:
#     flag_continue = False
#     ji_rho_in = get_indices_in_polygon(ilon_vec, ilat_vec, NR, NC)
#     dap = da[ji_rho_in[:,0], ji_rho_in[:,1]]
#     dvp = h[ji_rho_in[:,0], ji_rho_in[:,1]] * da[ji_rho_in[:,0], ji_rho_in[:,1]]
#     ap = np.nansum(dap)
#     vp = np.nansum(dvp)
#     hp = vp/ap
#     print('- all values are for unmasked area -')
#     print('Volume inside polygon = %0.1f km3' % (vp/1e9) )
#     print('Area inside polygon = %0.1f km2' % (ap/1e6) )
#     print('Mean Depth inside polygon = %0.1f m' % (hp) )
#     inp = input('Push Return to continue\n')
#     remove_poly(pline)
#     ax1.set_title('PAUSED')
# elif (bdict[nb]=='lineSave') and not flag_start:
#     flag_continue = False
#     pname = input('Name for saved polygon or line: ')
#     poutdir = Ldir['data'] / 'section_lines'
#     Lfun.make_dir(poutdir)
#     lon_poly = lonvec[ilon_vec]
#     lat_poly = latvec[ilat_vec]
#     pdict = {'lon_poly': lon_poly, 'lat_poly': lat_poly}
#     poutfn = poutdir / (pname.replace(' ','') + '.p')
#     pickle.dump(pdict, open(poutfn, 'wb'))
#     print(' - saved to ' + str(poutfn))
#     #remove_poly()
#     ax1.set_title('PAUSED')
#     if extra_map == True:
#         ax3.plot(lon_poly, lat_poly, '-*k', linewidth=1)
#         ax3.text(np.array(lon_poly).mean(),np.array(lat_poly).mean(),
#             pname, ha='center',va='center')
# elif (bdict[nb]=='lineInfo') and not flag_start:
#     flag_continue = False
#     x = ilon_vec
#     y = ilat_vec
#     dist = 0
#     for ii in range(len(x)-1):
#         dist += np.sqrt( (xm[y[ii+1],x[ii+1]]-xm[y[ii],x[ii]])**2
#             + (ym[y[ii+1],x[ii+1]]-ym[y[ii],x[ii]])**2 )
#     print('Line length = %0.1f km' % (dist/1e3))
#     inp = input('Push Return to continue\n')
#     remove_poly(pline)
#     ax1.set_title('PAUSED')


