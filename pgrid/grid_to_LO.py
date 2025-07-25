"""
This is the last step in the process, where you create and
populate LO_data/grids/[gridname].

This also makes the nudging to climatology NetCDF file.

"""
import numpy as np
import pickle
import shutil
import xarray as xr

from lo_tools import Lfun

import gfun
import gfun_utility as gfu

testing = True
if testing:
    from importlib import reload
    reload(gfu)

Gr =gfun.gstart()
# select grid file
in_fn = gfun.select_file(Gr)

# import gfun_user, to get s_dict
Ldir = Lfun.Lstart()
pth = Ldir['LO'] / 'pgrid'
upth = Ldir['LOu'] / 'pgrid'
if (upth / 'gfun_user.py').is_file():
    print('Importing gfun_user from LO_user')
    gfun_user = Lfun.module_from_file('gfun_user', upth / 'gfun_user.py')
else:
    print('Importing gfun_user from LO')
    gfun_user = Lfun.module_from_file('gfun_user', pth / 'gfun_user.py')

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# get grid size
ds = xr.open_dataset(in_fn)
NR, NC = ds.lon_rho.shape
ds.close()

N = gfun_user.s_dict['N']

# make output directory
Ldir = Lfun.Lstart(gridname=Gr['gridname'])
out_dir = Ldir['grid']
Lfun.make_dir(out_dir, clean=True)

# write s-coordinate info
s_out_fn = out_dir / 'S_COORDINATE_INFO.csv'
with open(s_out_fn, 'x') as f:
    f.write('ITEMS,VALUES\n')
    for k in gfun_user.s_dict.keys():
        f.write(k + ',' + str(gfun_user.s_dict[k]) + '\n')

# write xy-coordinate info
xy_out_fn = out_dir / 'XY_COORDINATE_INFO.csv'
with open(xy_out_fn, 'x') as f:
    f.write('NROWS,' + str(NR) + '\n')
    f.write('NCOLS,' + str(NC) + '\n')

# write dch info
dch_out_fn = out_dir / 'dch.csv'
with open(dch_out_fn, 'x') as f:
    for k in dch.keys():
        f.write(k + ',' + str(dch[k]) + '\n')

# copy files

shutil.copyfile(in_fn, out_dir / 'grid.nc')

try:
    shutil.copyfile(Gr['gdir'] / 'roms_river_info.csv', out_dir / 'river_info.csv')
except:
    # sometimes there is no river_info file
    pass

# Make the nudging to climatology file
gfu.make_nudgcoef(dch, out_dir, N, NR, NC)

# Run traps_placement
if 'do_traps' in dch.keys():
    if dch['do_traps']:
        print('\nRunning traps_placement and excluding selected items')
        from subprocess import Popen as Po
        from subprocess import PIPE as Pi
        import pandas as pd
        cmd_list = ['python', str(Ldir['LO'] / 'pre' / 'trapsP00' / 'traps_placement.py'), '-g', Gr['gridname']]
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stderr) > 0:
            print(stderr.decode())
        # if len(stdout) > 0:
        #     print(stdout.decode())
        # and edit out excluded triv and wwtp
        tw_dir = Ldir['data'] / 'grids' / Gr['gridname']
        for tw in ['triv','wwtp']:
            in_fn = tw_dir / (tw + '_info.csv')
            if in_fn.is_file():
                in_df = pd.read_csv(in_fn, index_col='rname')
                for item in dch['excluded_' + tw]:
                    in_df = in_df.drop(item)
                    print('Excluding %s' % (item))
                in_df.to_csv(in_fn)

