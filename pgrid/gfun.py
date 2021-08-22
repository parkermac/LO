# -*- coding: utf-8 -*-
"""
Organizational functions for pgrid.
"""

# **** USER EDIT ********
#gridname = 'aestus3'
gridname = 'cas6'
# **** END USER EDIT ****

import os; import sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))

dir0 = Ldir['parent']
pgdir = dir0 + 'ptools_output/pgrid/'

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
    dch['t_dir'] = Gr['dir0'] + 'ptools_data/topo/'    
    # list of topo files: coarsest to finest
    dch['t_list'] = ['srtm15/topo15.nc',
              'cascadia/cascadia_gridded.nc',
             'psdem/PS_183m.nc',
             'ttp_patch/TTP_Regional_27m_patch.nc']
 
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

def gstart(gridname=gridname):
    
    if gridname in ['aestus1', 'aestus2']:
        ri_dir = dir0 + 'ptools_output/river/analytical/'
    else:
        ri_dir = dir0 + 'ptools_output/river/pnw_all_2016_07/'
    
    gdir = pgdir + gridname + '/'
    Gr ={'gridname': gridname, 'dir0': dir0, 'pgdir': pgdir, 'gdir': gdir,
         'ri_dir': ri_dir}
    return Gr

def select_file(Gr, using_old_grid=False):
    # interactive selection
    if using_old_grid==True:
        fn_list = []
        dir0 = Ldir['parent'] + 'LiveOcean_data/grids/'
        gn_list = ['cascadia1', 'cascadia2']
        for gn in gn_list:
            fn_list.append(dir0 + gn + '/grid.nc')
    elif using_old_grid==False:
        print('\n** %s in <<%s>> **\n' % ('Choose file to edit', Gr['gridname']))
        fn_list_raw = os.listdir(Gr['gdir'])
        fn_list = []
        for item in fn_list_raw:
            if item[-3:] == '.nc':
                fn_list.append(item)
        fn_list.sort()
    Nfn = len(fn_list)
    fn_dict = dict(zip(range(Nfn), fn_list))
    for nfn in range(Nfn):
        print(str(nfn) + ': ' + fn_list[nfn])
    my_nfn = int(input('-- Input number -- '))
    fn = fn_dict[my_nfn]
    return fn

def increment_filename(fn, tag='_m'):
    # create the new file name
    gni = fn.find(tag)
    new_num = ('00' + str(int(fn[gni+2: gni+4]) + 1))[-2:]
    fn_new = fn.replace(fn[gni:gni+4], tag + new_num)

    return fn_new
