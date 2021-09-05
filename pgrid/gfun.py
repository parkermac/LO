"""
Organizational functions for pgrid.
"""

# **** USER EDIT ********

# This is the name of the grid that you are working on.
gridname = 'sal0'

# These are the gridname and tag to feed to use when creating the Ldir paths.
# They are used for accessing the river tracks, which may be developed for one
# grid but reused in others.
base_gridname = 'cas6'
base_tag = 'v3'

# **** END USER EDIT ****

from pathlib import Path
from lo_tools import Lfun

Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)

def gstart():
    """
    This returns a dict of Path objects that tell where various things are,
    or where they should go.
    """
    pgdir = Ldir['LOo'] / 'pgrid'
    gdir = pgdir / gridname # where grid.nc will end up
    ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'tracks'
    Gr ={'gridname': gridname,'pgdir': pgdir, 'gdir': gdir,'ri_dir': ri_dir}
    return Gr

def default_choices(Gr, wet_dry=False):
    # Default choices (can override in each case)    
    dch = dict()    

    # Decide if the grid will allow wetting and drying.
    # We do this first becasue it affects several subsequent choices
    dch['wet_dry'] = wet_dry

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
    if dch['wet_dry'] == True:
        dch['unmask_coast'] = True
    else:
        dch['unmask_coast'] = False
    # Set remove_islands to True to automatically remove isolated patches of
    # land or ocean.
    dch['remove_islands'] = True

    # SMOOTHING
    dch['use_min_depth'] = True # now I think this is always a good idea
    dch['min_depth'] = 4 # meters (positive down)
        
    # NUDGING
    # Use nudging edges to decide which edges to have nudging to climatology
    # on. And the nudging_days are the the (short, long) timescales to use.
    dch['nudging_edges'] = ['north', 'south', 'east', 'west']
    dch['nudging_days'] = (3.0, 60.0)
        
    return dch

def select_file(Gr):
    # interactive selection
    fn = Lfun.choose_item(Gr['gdir'], tag='.nc', itext='** Choose grid from list **')
    return fn

def increment_filename(fn, tag='_m'):
    fns = str(fn)
    # create the new file name
    gni = fns.find(tag)
    new_num = ('00' + str(int(fns[gni+2: gni+4]) + 1))[-2:]
    fns_new = fns.replace(fns[gni:gni+4], tag + new_num)
    fn_new = Path(fns_new)
    return fn_new
