#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 09:15:30 2017

@author: PM5

This is a convenience program to allow you to edit the items in the
dch "choices" dictionary, without having to start a new version of the grid.

This should be used sparingly becasue it circumvents the logic of making
all imortant choices up-front in make_grid.py.  However it can be useful
for jobs like changing the nudging times used by make nudgcoef.
"""

import pickle
import matplotlib.pyplot as plt
import netCDF4 as nc

from importlib import reload
import gfun; reload(gfun)
import gfun_utility as gfu
reload(gfu)
Gr =gfun.gstart()
import Lfun
Ldir = Lfun.Lstart(gridname=Gr['gridname'])
out_dir = Ldir['grid']

def plot_nudgcoef(gridname, out_dir):
    ds = nc.Dataset(out_dir + 'nudgcoef.nc')    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap(name='jet')
    nudc = ds['M2_NudgeCoef'][:]
    cs = ax.pcolormesh(nudc, cmap = cmap)
    fig.colorbar(cs, ax=ax, extend='both')
    ax.set_title(gridname + ' Nudging timescales (1/day)')
    ds.close()
    plt.show()
    return nudc
    
# load the default choices
dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

# update the choices

if Gr['gridname'] == 'sal0':
    # changing nudging parameters
    dch['nudging_edges'] = ['north', 'west']
    dch['nudging_days'] = (0.1, 1.5)
    # and add the revised file to LiveOcean_data/grids
    # without cleaning out the directory
    gfu.make_nudgcoef(dch, out_dir)
    nudc = plot_nudgcoef(Gr['gridname'], out_dir)
    L2 = int(nudc.shape[0]/2)
    M2 = int(nudc.shape[1]/2)
    print('North: ' + 10*'%0.1f, ' % tuple([i for i in nudc[-10:,L2]]))
    print('South: ' + 10*'%0.1f, ' % tuple([i for i in nudc[:10,L2]]))
    print('East: ' + 10*'%0.1f, ' % tuple([i for i in nudc[M2,-10:]]))
    print('West: ' + 10*'%0.1f, ' % tuple([i for i in nudc[M2,:10]]))


# save the default choices
pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))

# also save the dch dict
dch_fn = out_dir + 'dch.csv'
Lfun.dict_to_csv(dch, dch_fn)


