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
start_line_ax = plt.subplot2grid((4,4), (0,3))
done_ax = plt.subplot2grid((4,4), (1,3))
tb_name_ax = plt.subplot2grid((4,4), (2,3))
tb_remove_ax = plt.subplot2grid((4,4), (3,3))
# the widgets
start_line_button = Button(start_line_ax, 'Start Line', color='0.85', hovercolor='0.95')
done_button = Button(done_ax, 'Done', color='0.85', hovercolor='0.95')
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

W = dict()
W['done'] = False

def name_section(sect_name):
    ax1.set_title('Click to Add Points (Return to Stop)')
    plt.draw()
    df = pd.read_pickle(temp_dir / 'temp_df.p')
    ilon_vec = df.ilon.to_numpy()
    ilat_vec = df.ilat.to_numpy()
    if len(sect_name) > 0:
        df.to_pickle(temp_dir / (sect_name + '.p'))
    ax1.text(ilon_vec[0], ilat_vec[0], sect_name)
    tb_name.set_val('')
    ax1.set_title('PAUSED', color='g')
    plt.draw()

def on_start_line(junk):
    ax1.set_title('Click to Add Points (Return to Stop)')
start_line_button.on_clicked(on_start_line)

def on_done(junk):
    W['done'] = True
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
    sys.exit()
done_button.on_clicked(on_done)

def remove_section(x_sect_name):
    dfr = pd.read_pickle(c_fn)
    xx = dfr.sn==x_sect_name
    dfr = dfr.drop(index=dfr.index[xx], axis=1)
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

plt.show()


while W['done'] == False:
    a = plt.ginput(n=10, timeout=0)
    print(len(a))
    # a is a list of tuples
    ilon_list = []
    ilat_list = []
    for b in a:
        ix = np.round(b[0]).astype(int)
        iy = np.round(b[1]).astype(int)
        ilon_list.append(ix)
        ilat_list.append(iy)
    ilon_vec = np.array(ilon_list)
    ilat_vec = np.array(ilat_list)
    aa = ax1.axis()
    ax1.plot(ilon_vec, ilat_vec,'*-r')
    ax1.axis(aa)
    ax1.set_title('Done with that line')
    plt.draw()

    temp_df = pd.DataFrame(columns=['ilon','ilat'])
    temp_df.loc[:,'ilon'] = ilon_vec
    temp_df.loc[:,'ilat'] = ilat_vec
    temp_df.to_pickle(temp_dir / 'temp_df.p')
    ax1.set_title('Saved to temp_df.p')
    tb_name.on_submit(name_section)
    



        