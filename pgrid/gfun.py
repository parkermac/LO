"""
Organizational functions for pgrid.
"""

from lo_tools import Lfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO_user' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_user

# assume the user will be editing things
from importlib import reload
reload(gfun_user)

# get info from LO_user/pgrid/gfun_user.py
gridname = gfun_user.gridname
base_gridname = gfun_user.base_gridname
base_tag = gfun_user.base_tag

Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)

def gstart(gridname=gridname):
    """
    This returns a dict of Path objects that tell where various things are,
    or where they should go.
    """
    pgdir = Ldir['LOo'] / 'pgrid'
    gdir = pgdir / gridname # where grid.nc will end up
    ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
    Gr ={'gridname': gridname,'pgdir': pgdir, 'gdir': gdir,'ri_dir': ri_dir}
    return Gr

def default_choices():
    # Default choices (can override in each case)
    dch = dict()

    # Decide if the grid will allow wetting and drying.
    # We do this first because it affects several subsequent choices
    # dch['wet_dry'] = wet_dry
    # deprecated 2021.11.15

    # GRID CREATION
    # Set analytical to true when we define the bathymetry analytically.
    dch['analytical'] = False
    
    # z_offset is an adjustment to zero of the bathymetry to account for
    # the fact that mean sea level is somewhat higher than NAVD88.
    dch['use_z_offset'] = True
    dch['z_offset'] = -1.06
    
    # specify topography files to use
    t_dir = Ldir['data'] / 'topo'
    dch['t_dir'] = Ldir['data'] / 'topo'
    # list of topo files: coarsest to finest
    dch['t_list'] = [t_dir / 'srtm15' / 'topo15.nc',
              t_dir / 'cascadia' / 'cascadia_gridded.nc',
             t_dir / 'psdem' / 'PS_183m.nc',
             t_dir / 'ttp_patch' / 'TTP_Regional_27m_patch.nc']
 
    # MASKING
    # list of existing masks to work from
    dch['maskfiles'] = []
    # set z position of INITIAL dividing line (positive up)
    dch['z_land'] = 0
    # Set unmask_coast to True to unmask all cells crossed by the coastline.
    dch['unmask_coast'] = False
    # Set remove_islands to True to automatically remove isolated patches of
    # land or ocean.
    dch['remove_islands'] = True

    # SMOOTHING
    dch['use_min_depth'] = True
    dch['min_depth'] = 4 # meters (positive down)
        
    # NUDGING
    # Use nudging edges to decide which edges to have nudging to climatology
    # on. And the nudging_days are the the (short, long) timescales to use.
    dch['nudging_edges'] = ['north', 'south', 'east', 'west']
    dch['nudging_days'] = (3.0, 60.0)
        
    return dch

def select_file(Gr):
    # interactive selection
    fn = Lfun.choose_item(Gr['gdir'], tag='.nc',
        itext='** Choose grid from list **', last=True)
    return  Gr['gdir'] / fn

def increment_filename(fn, tag):
    """
    Creates an updated filename (Path object), increasing the number
    for the specified tag.
    """
    if tag not in ['_m', '_r', '_s', '_x']:
        print('Error in increment_filename: unrecognized tag.')
        return
    fns = fn.name
    # create the new file name
    gni = fns.find(tag)
    new_num = ('00' + str(int(fns[gni+2: gni+4]) + 1))[-2:]
    fns_new = fns.replace(fns[gni:gni+4], tag + new_num)
    fn_new = fn.parent / fns_new
    return fn_new
