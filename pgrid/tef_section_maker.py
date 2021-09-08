# -*- coding: utf-8 -*-
"""
Tool for interactively making TEF sections.  The screen output can be copied
as lines to be used in tef_fun.py

Can be run in ipython with a user-specified grid file

run tef_section_maker.py -g sal0

"""

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
import pickle

import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', default='',
        type=str)
args = parser.parse_args()
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
    
# select grid file
using_old_grid = False
# Set this to True to look at grids we have already created,
# e.g. ones currently in use for LiveOcean.
# Set it to False when interacting with grids from pgrid_output.
if using_old_grid==True:
    fn = gfun.select_file(Gr, using_old_grid=True)
    in_fn = fn
elif using_old_grid==False:
    fn = gfun.select_file(Gr)
    in_fn = Gr['gdir'] + fn

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
fig = plt.figure(figsize=(22,12)) # (13,8) is good for my laptop
ax1 = plt.subplot2grid((1,5), (0,0), colspan=2) # imshow map
ax2 = plt.subplot2grid((1,5), (0,2), colspan=1) # buttons
ax3 = plt.subplot2grid((1,5), (0,3), colspan=2) # lat-lon map

#%% initialize the data plot
cmap1 = plt.get_cmap(name='rainbow_r') # terrain
tvmin = -20
tvmax = 200
cs = ax1.imshow(h, interpolation='nearest', vmin=tvmin, vmax=tvmax, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
aa = ax1.axis()
# add the coastline
clon, clat = pfun.get_coast()
cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:])
cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0])
ax1.plot(cx0 + cxf, NR - (cy0 + cyf) - 1, '-k')
# add rivers
gfp.edit_mask_river_tracks(Gr, NR, ax1)
ax1.axis(aa)

# create control buttons
# list is organized from bottom to top
blist = ['start', 'pause', 'startSect',
         'sectSave', 'done']
# nicer names
Blist = ['Start', 'Pause', 'Start Section',
         'Save Section', 'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

# initialize the lat,lon map plot
ax3.pcolormesh(plon, plat, Hm,
    vmin=tvmin, vmax=tvmax, cmap = cmap1)
pfun.add_coast(ax3)
ax3.axis(aa0)
pfun.dar(ax3)
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')

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
        elif (bdict[nb]=='startSect') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Begin Section')
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='sectSave') and not flag_start:
            flag_continue = False
            pname = input('Name for saved section: ')
            landward = input('Sign for landward (1 or -1): ')
            lon_poly = lonvec[plon_poly]
            lat_poly = latvec[plat_poly]
            ax1.set_title('PAUSED')
            ax3.plot(lon_poly, lat_poly, '-*k', linewidth=1)
            ax3.text(np.array(lon_poly).mean(),np.array(lat_poly).mean(),pname)
            # print output for use in tef_fun.py
            tef_list_item = ('sect_df.loc[\'%s\',:] = [%8.3f, %8.3f, %8.3f, %8.3f, %d]' %
                (pname, lon_poly[0], lon_poly[1], lat_poly[0], lat_poly[1], int(landward)) )
            tef_list.append(tef_list_item)
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE', fontweight='bold', color='b')
            
            print('\n*** Entries for tef_fun.get_sect_df() **\n')
            if len(tef_list) > 0:
                for item in tef_list:
                    print(item)
        else:
            pass
        plt.draw()

    # this code deals with map input, and only responds when
    # we are clicking on the map
    elif flag_continue and not flag_start:
        # we are in the data field
        if (flag_e == 'p') and (len(plon_poly) < 2):
            # this draws a polygon as you click
            plon_poly.append(ix)
            plat_poly.append(iy)
            remove_poly()
            aa = ax1.axis()
            pline = ax1.plot(plon_poly, plat_poly,'*-r')
            ax1.axis(aa)
            ax1.set_title('Add to Section')
            if (len(plon_poly) == 2):
                flag_continue = False
                ax1.set_title('PAUSED')
                # adjust indices to make it perfectly zonal or meridional
                dix = np.abs(plon_poly[1] - plon_poly[0])
                diy = np.abs(plat_poly[1] - plat_poly[0])
                if dix > diy: # EW section
                    plat_poly[1] = plat_poly[0]
                elif diy > dix: # NS section
                    plon_poly[1] = plon_poly[0]
        plt.draw()


