"""
This is the last step in the process, where you create and
populate LO_data/grids/[gridname].

This also makes the nudging to climatology NetCDF file.

"""
import numpy as np
import pickle
import shutil
# import xarray as xr

from lo_tools import Lfun

import gfun
import gfun_utility as gfu

Gr =gfun.gstart()
# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_x')

# # get the grid from NetCDF
# ds = xr.open_dataset(in_fn)
# h = ds.h.values
# mask_rho = ds.mask_rho.values

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))
# ---

# copy to LiveOcean_data/grids

Ldir = Lfun.Lstart(gridname=Gr['gridname'])

# make output directory
out_dir = Ldir['grid']
Lfun.make_dir(out_dir, clean=True)

# copy files
shutil.copyfile(in_fn, out_dir / 'grid.nc')
try:
    shutil.copyfile(Gr['gdir'] / 'roms_river_info.csv', out_dir / 'river_info.csv')
except:
    # sometimes there is no river_info file
    pass

# also save the dch dict
dch_fn = out_dir / 'dch.csv'
Lfun.dict_to_csv(dch, dch_fn)

# #%% and put in the S coordinate info
#
# # copy file
# if False:
#     # regular 40-layer grid
#     shutil.copyfile(Ldir['data'] + 'grids/S_COORDINATE_INFO_1.csv',
#                     out_dir + 'S_COORDINATE_INFO.csv')
# else:
#     # new 30-layer grid
#     print('WARNING: using 30 layer grid!!')
#     shutil.copyfile(Ldir['data'] + 'grids/S_COORDINATE_INFO_3.csv',
#                     out_dir + 'S_COORDINATE_INFO.csv')
#
# # create the S.mat file
# import subprocess
# func = ( "Z_make_S(\'" + out_dir + "\')" )
# cmd = Ldir['which_matlab']
# run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
# # I dropped the "-nojvm" because it was throwing an error.  Maybe
# # this is because I am using a newer matlab release.
# cwd = Ldir['LO']+'shared/Z_functions/'
# # I had to add the cwd= parameter in order to run in that directory.
# proc = subprocess.Popen(run_cmd, cwd=cwd,
#                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# out, err = proc.communicate()
# # "out" is the screen output of the matlab code, and
# # "err" stores arror messages from subprocess (I think).
# #print(out.decode()) # will print "out" nicely

# Make the nudging to climatology file
gfu.make_nudgcoef(dch, out_dir)
