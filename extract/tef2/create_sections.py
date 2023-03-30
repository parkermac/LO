"""
Code to interactively create a collection of TEF sections for
a user-specified grid. I am trying to keep it simple, using
a limited number of keyboard inputs for all actions.

Command line arguments:

-g, --gridname: which grid.nc to use
-ctag, --collection_tag: save collection in the folder
    LO_output/extract/tef2/sections_[gctag]
    where [gctag] = [gridname]_[collection_tag]
-clobber: True to start clean collection folder
-small: True to work on a laptop

Examples:

run create_sections -g cas6
This would start or append to a [gctag] called cas6_test.

You can look at an existing [gctag] with (for me) a command like:
run create_sections -g cas6 -ctag c0
and then just do "q" to exit without changing anything

"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cmocean import cm
import sys

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# gridname and tag will create a collection folder
parser.add_argument('-g','--gridname', default='ae0', type=str)
parser.add_argument('-ctag','--collection_tag', default='test', type=str)
# set clobber to True to start a new collection even if it exists
parser.add_argument('-clobber', default=False, type=Lfun.boolean_string)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input and output locations
Ldir = Lfun.Lstart(gridname=args.gridname)
in_fn = Ldir['grid'] / 'grid.nc'
out_name = 'sections_' + args.gridname + '_' + args.collection_tag
out_dir = Ldir['LOo'] / 'extract' / 'tef2' / out_name
if args.clobber:
    inp = input('Do you really want to clobber? y/n: ')
    if inp == 'y':
        Lfun.make_dir(out_dir, clean=True)
    else:
        sys.exit()
else:
    # make out_dir in case it does not exist yet
    Lfun.make_dir(out_dir)
    
# get grid data
ds = xr.open_dataset(in_fn)
h = ds.h.values
m = ds.mask_rho.values
h[m==0] = np.nan
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
ds.close

# initialize plot
plt.close('all')
if args.small == True:
    fig = plt.figure(figsize=(8,8)) # laptop size
else:
    fig = plt.figure(figsize=(12,12)) # external monitor size
ax = fig.add_subplot(111)
ax.pcolormesh(plon,plat,h, vmin=-30, vmax=200, cmap=cm.deep)
pfun.dar(ax)
ax.set_title('Use keyboard to create or remove sections')
ax.text(.05,.95,out_name,transform=ax.transAxes,
    fontweight='bold',bbox=pfun.bbox)
plt.show()

# dicts to keep line and text handles
ld = dict(); td = dict()

def plot_line(sn):
    df = pd.read_pickle(out_dir / (sn + '.p'))
    x = df.x.to_numpy()
    y = df.y.to_numpy()
    ld[sn] = ax.plot(x,y,'-or', lw=2)
    td[sn] = ax.text(x[0],y[0],'\n'+sn,c='r',ha='center',va='top',
        fontweight='bold')
    plt.draw()

def get_sn_list():
    df_list = list(out_dir.glob('*.p'))
    sn_list = [item.name.replace('.p','') for item in df_list]
    return sn_list

if args.clobber == False:
    # plot all existing sections in the collection
    sn_list = get_sn_list()
    for sn in sn_list:
        plot_line(sn)

# interactive section creation
print('=========== tef2 section creation tool =========================')
print('You need to go back and forth between the plot and the keyboard.')
print('When creating a section you need to click on the plot to make it active.')
print('The cursor will then appear as a cross.')
print('It is OK to resize anytime while waiting for keyboard input.')
print('-> Example command to create a section called a1: c a1')
print('-> Example command to delete a section called a1: x a1')
print('-> Type q to quit')
print('Sections are saved as pickled pandas DataFrames in:')
print(out_dir)
print('')
keep_going = True
while keep_going:
    
    # parse keyboard input
    inp = input('Input choices c [],x [],q: ')
    ii = inp.split(' ')
    if len(ii) == 1:
        if ii[0] != 'q':
            print('Error: q is the only single-character input, try again.')
            task = 'junk'
        else:
            task = 'q'
    elif len(ii) == 2:
        task = ii[0]
        sn = ii[1] # section name
        sn_list = get_sn_list()
        if (sn in sn_list) and (task == 'c'):
            print('Oops, that name is already used. Please delete it first.')
            task = 'junk'
    else:
        print('Error reading keyboard input, try again (q to quit).')
        task = 'junk'
    
    # do the requested task
    if task == 'c':
        ax.set_title('Click to create section, Return to end')
        plt.draw()
        p = fig.ginput(n=-1) # returns a list of tuples
        x = []; y = []
        x = [item[0] for item in p]
        y = [item[1] for item in p]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = x
        df.loc[:,'y'] = y
        df.to_pickle(out_dir / (sn + '.p'))
        plot_line(sn)
    elif task == 'x':
        try:
            ld[sn][0].set_visible(False)
            td[sn].remove()
        except KeyError:
            print(sn + ' not in dicts of lines or labels')
        sn_list = get_sn_list()
        if sn in sn_list:
            (out_dir / (sn + '.p')).unlink(missing_ok=True)
        else:
            print(sn + ' not among pickled lines')
        plt.draw()
    elif task == 'q':
        keep_going = False
        ax.set_title('DONE')
        plt.draw()
    else:
        pass
