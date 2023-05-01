"""
Code to interactively create a collection of river tracks for
a user-specified grid. I am trying to keep it simple, using
a limited number of keyboard inputs for all actions.

This is specialized code in case you want to create new river tracks
or change ones you alreay have.

NOTES:

(1) The standard treatment of these in subsequent river code assumes
that they start from the ocean end.

(2) You have to already have a grid to start from, which is a little
clumsy, but will do for now. You could just make one with the right
region as a start. You have to go all the way through the final step
"grid_to_LO" so that there will be a grid file to find. The tracks you
create here will not change anything until you run carve_rivers.

(3) You still have to do some work to make these function right in
the carve_rivers step, and in subsequent river forcing steps.

Command line arguments:

-g, --gridname: which grid.nc to use
-ctag, --collection_tag: save collection in the folder
    LO_output/pre/river1/[ctag]
-clobber: True to start clean collection folder
-small: True to work on a laptop

Examples:

run create_river_tracks -g wgh0 -ctag lo_new
This would start or append to a [ctag] called lo_new.

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
parser.add_argument('-g','--gridname', default='wgh0', type=str)
parser.add_argument('-ctag','--collection_tag', default='lo_new', type=str)
# set clobber to True to start a new collection even if it exists
parser.add_argument('-clobber', default=False, type=Lfun.boolean_string)
# set small to True to work on a laptop
parser.add_argument('-small', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# input and output locations
Ldir = Lfun.Lstart(gridname=args.gridname)
in_fn = Ldir['grid'] / 'grid.nc'
ctag = args.collection_tag
out_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag / 'tracks'
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
ax.pcolormesh(plon,plat,h, vmin=0, vmax=15, cmap=cm.deep)
aa = ax.get_xlim() + ax.get_ylim()
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_title('Use keyboard to create or remove tracks')
ax.text(.05,.95,ctag,transform=ax.transAxes,
    fontweight='bold',bbox=pfun.bbox)
plt.show()

# dicts to keep line and text handles
ld = dict(); td = dict()

def plot_line(sn):
    df = pd.read_pickle(out_dir / (sn + '.p'))
    x = df.lon.to_numpy()
    y = df.lat.to_numpy()
    ld[sn] = ax.plot(x,y,'-or', lw=2)
    td[sn] = ax.text(x[-1],y[-1],'\n'+sn,c='r',ha='left',va='bottom',
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
print('=========== river track creation tool =========================')
print('You need to go back and forth between the plot and the keyboard.')
print('When creating a track you need to click on the plot to make it active.')
print('The cursor will then appear as a cross.')
print('It is OK to resize anytime while waiting for keyboard input.')
print('-> Example command to create a track called willapa: c willapa')
print('-> Example command to delete a track called willapa: x willapa')
print('-> Type q to quit')
print('Tracks are saved as pickled pandas DataFrames in:')
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
        ax.set_title('Click to create a river track, starting from the ocean, Return to end')
        plt.draw()
        p = fig.ginput(n=-1, timeout=-1, mouse_stop=None) # returns a list of tuples
        # Note that using mouse_stop=None in the command above solves a problem
        # I was having where a stray touch of the mouse would cause termination of the
        # line. Now you can only terminate by hitting return.
        x = []; y = []
        x = [item[0] for item in p]
        y = [item[1] for item in p]
        df = pd.DataFrame(columns=['lon','lat'])
        df.loc[:,'lon'] = x
        df.loc[:,'lat'] = y
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
