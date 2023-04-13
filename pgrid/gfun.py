"""
Organizational functions for pgrid.
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# get the gfun_user module, looking first in LO_user
pth = Ldir['LO'] / 'pgrid'
upth = Ldir['LOu'] / 'pgrid'
if (upth / 'gfun_user.py').is_file():
    print('Importing gfun_user from LO_user')
    gfun_user = Lfun.module_from_file('gfun_user', upth / 'gfun_user.py')
else:
    print('Importing gfun_user from LO')
    gfun_user = Lfun.module_from_file('gfun_user', pth / 'gfun_user.py')

# get info from LO_user/pgrid/gfun_user.py
gridname = gfun_user.gridname

# remake Ldir
Ldir = Lfun.Lstart()

def gstart(gridname=gridname):
    """
    This returns a dict of Path objects that tell where various things are,
    or where they should go.
    """
    pgdir = Ldir['LOo'] / 'pgrid'
    gdir = pgdir / gridname # where grid.nc will end up
    ri_dir0 = Ldir['LOo'] / 'pre' / 'river1'
    Gr ={'gridname': gridname,'pgdir': pgdir, 'gdir': gdir, 'ri_dir0': ri_dir0}
    return Gr

def default_choices():
    # Default choices (can override in each case)
    dch = dict()

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
    
    # RIVERS
    # ctag for river info, using the new LO/pre/river1 system
    dch['ctag'] = 'lo_base'
    # Sometimes for nests we mask out a partial basin, and this allows us to also
    # exclude its rivers.  The downside is that you need to know what those rivers
    # are in advance.  If you want to do this you should run carve_rivers.py and
    # decide what to exclude, and then start again, before putting a lot of
    # time in to edit_mask.py.
    dch['excluded_rivers'] = []

    # SMOOTHING
    dch['use_min_depth'] = True
    dch['min_depth'] = 4 # meters (positive down)
    dch['rx0max'] = 0.15
        
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
