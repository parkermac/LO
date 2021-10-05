"""
Code for particle tracking, designed for ROMS output.  This new version
makes extensive use of nearest-neighbor KDTree algorithms for interpolation.
This results is significantly (36x) faster runtimes compared with old version.

PERFORMANCE: about 3 minutes per day for a 3D cas6 experiment with 10k particles.

NOTE: You have to have run make_KDTrees.py for the grid (e.g. cas6) before running.

NOTE: There is some issue, perhaps with garbage collection, which causes
the loading of NetCDF files to happen slower after running a few times
interactively from ipython.  It appears that this can be avoided by running
from the terminal as: python tracker.py [args].

This program is a driver where you specify:
- an experiment (ROMS run + release locations + other choices)
- a release or set of releases within that experiment (start day, etc.)

The main argument you provide is -exp, which is the experiment name, and
is used by experiments.get_exp_info() and .get_ic() to get the gtagex and initial particle
locations.  Other possible commmand line arguments and their defaults
are explained in the argparse section below.

NOTE: To improve usefulness for people other than me, this driver will
first look for:
- LiveOcean_user/tracker/user_trackfun.py
before loading my versions.

This allows you to create your own modifications to the tracking
(e.g. for diurnal depth behavior) while still being able to use git pull to update the main code.

It can be run on its own, or with command line arguments to facilitate
large, automated jobs, for example in python:
    
[examples under construction]

"""

import sys
from datetime import datetime, timedelta
from time import time
import argparse
import numpy as np

from lo_tools import Lfun, zfun
Ldir = Lfun.Lstart()

from importlib import reload

pth = Ldir['LOu'] / 'tracker'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import experiments as exp
reload(exp)

import trackfun_nc as tfnc
reload(tfnc)

# The import of trackfun or user_trackfun is done later in this program,
# about 100 lines down.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
# command line arguments, can be input in any order
parser = argparse.ArgumentParser()

# Set the experiment name
# (details set in experiments.py, or, if it exists, user_experiments.py)
parser.add_argument('-exp', '--exp_name', default='jdf0', type=str)
parser.add_argument('-clb', '--clobber', default=False, type=zfun.boolean_string)
# overwrite existing output folder if clobber == True
parser.add_argument('-t', '--tag', default='', type=str)
# append an optional tag to the end of the output folder name

# These are False unless the flags are used with the argument True
# so if you do NOT use these flags the run will be:
# - trapped to the surface
# - no vertical turbulent diffusion
parser.add_argument('-3d', default=False, type=zfun.boolean_string) # do 3d tracking
parser.add_argument('-laminar', default=False, type=zfun.boolean_string) # no turbulence
parser.add_argument('-no_advection', default=False, type=zfun.boolean_string) # no advection
parser.add_argument('-sink', default=0, type=float) # particle sinking speed (m per day, e.g. 40)

# windage = a small number: 0 <= windage << 1 (e.g. 0.03)
# fraction of windspeed added to advection, only for 3d=False
parser.add_argument('-wnd', '--windage', default=0, type=float)

# set the starting day
parser.add_argument('-ds', '--ds_first_day', default='2019.07.04', type=str)

# You can make multiple releases using:
# number_of_start_days > 1 & days_between_starts, and which hour (UTC) to start on
parser.add_argument('-nsd', '--number_of_start_days', default=1, type=int)
parser.add_argument('-dbs', '--days_between_starts', default=1, type=int)
parser.add_argument('-dtt', '--days_to_track', default=1, type=int)
parser.add_argument('-sh', '--start_hour', default=0, type=int)

# number of divisions to make between saves for the integration
# e.g. if ndiv = 12 and we have hourly saves, we use a 300 sec step
# for the integration. 300 s seems like a good default value,
# based on Banas et al. (2009, CSR RISE paper).
parser.add_argument('-ndiv', default=12, type=int)
parser.add_argument('-sph', default=1, type=int)
# sph = saves per hour, a new argument to allow more frequent writing of output.

# set which ROMS output directory to look in.  This allows you to choose
# between Ldir['roms'] or Ldir['roms2'] by passing "roms" or "roms2", or any other
# choice that eventually appears in Ldir.
parser.add_argument('-ro', '--roms_out_num', type=int, default=2) # 1 = Ldir['roms_out1'], etc.
# valid arguments to pass are: roms, roms2, roms3

args = parser.parse_args()
TR = args.__dict__ 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# set where to look for model output
if args.roms_out_num == 0:
    TR['roms_out'] = Ldir['roms_out']
elif args.roms_out_num > 0:
    TR['roms_out'] = Ldir['roms_out' + str(args.roms_out_num)]

# set dependent and default fields

TR['turb'] = False

# make sure sph is no greater than ndiv
TR['sph'] = np.min((TR['sph'],TR['ndiv']))

# overrides
if TR['3d']:
    TR['windage'] = 0
    TR['turb'] = True # default is that 3d is always turbulent
    
if TR['laminar']:
    TR['turb'] = False

# get experiment info
EI = exp.get_exp_info(TR['exp_name'])
TR['gtagex'] = EI['gtagex']
TR['gridname'] = EI['gridname']

# pass some info to Ldir
Ldir['gtagex'] = TR['gtagex']
Ldir['roms_out'] = TR['roms_out']

# get the full path to a valid history file
fn00 = Ldir['roms_out'] / TR['gtagex'] / ('f' + TR['ds_first_day']) / 'ocean_his_0001.nc'
TR['fn00'] = fn00
EI['fn00'] = str(fn00)

# set the name of the output folder
out_name = TR['exp_name']

# modify the output folder name, based on other choices
if TR['3d']:
    out_name += '_3d'
elif not TR['3d']:
    out_name += '_surf'
if TR['laminar']:
    out_name += '_laminar'
if TR['windage'] > 0:
    out_name += '_wind' + str(int(100*TR['windage']))
if TR['start_hour'] > 0:
    out_name += '_sh' + str(int(TR['start_hour']))
if TR['sink'] > 0:
    out_name += '_sink' + str(int(TR['sink']))
if TR['no_advection'] == True:
    out_name += '_nadv'
if TR['ndiv'] != 12: # only mention ndiv if it is NOT 12
    out_name += '_ndiv' + str(TR['ndiv'])
if len(TR['tag']) > 0:
    out_name += '_' + TR['tag']

# make the list of start days (datetimes) for separate releases
idt_list = []
dt = datetime.strptime(TR['ds_first_day'], '%Y.%m.%d')
for nic in range(TR['number_of_start_days']):
    idt_list.append(dt)
    dt = dt + timedelta(TR['days_between_starts'])

# make sure the output parent directory "tracks" exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = Ldir['LOo'] / 'tracks'
Lfun.make_dir(outdir0)

# make the output directory (empty)
outdir1 = out_name
outdir = outdir0 / outdir1
if outdir.is_dir():
    if args.clobber:
        pass # continue and overwrite if clobber is True
    else:
        print('Warning: output directory exists - rename if you want to keep it!!')
        print('-- tracker run not started --')
        sys.exit()
Lfun.make_dir(outdir, clean=True)
print(50*'*' + '\nWriting to ' + str(outdir))
sys.stdout.flush()

# and write some info to outdir0 for use by trackfun.py
ei_out_fn = outdir0 / 'exp_info.csv'
ei_out_fn.unlink(missing_ok=True)
with open(ei_out_fn, 'w') as f:
    for k in TR.keys():
        f.write(k + ',' + str(TR[k]) + '\n')
ei_out_fn = outdir / 'exp_info.csv'
ei_out_fn.unlink(missing_ok=True)
with open(ei_out_fn, 'w') as f:
    for k in TR.keys():
        f.write(k + ',' + str(TR[k]) + '\n')

# Load the trackfun module.
# NOTE: we have to load this module AFTER we write [outdir0]/exp_info.csv
# because it uses that information to decide which KDTrees to load.  Crude.
if (Ldir['LOu'] / 'tracker' / 'user_trackfun.py').is_file():
    sys.path.append(str(Ldir['LOu'] / 'tracker'))
    import user_trackfun as tfun
else:
    import trackfun as tfun
reload(tfun)

# get the initial particle location vectors
plon00, plat00, pcs00 = exp.get_ic(EI, fn00)

# step through the releases, one for each start day
write_grid = True
for idt0 in idt_list:
    tt0 = time() # monitor integration time
    
    # name the release file by start day
    idt0_str = datetime.strftime(idt0,'%Y.%m.%d')
    outname = ('release_' + idt0_str + '.nc')
    print('-- ' + outname)
    sys.stdout.flush()
    out_fn = outdir / outname
    
    # we do the calculation in one-day segments, but write complete
    # output for a release to a single NetCDF file.
    for nd in range(TR['days_to_track']):
        
        # get or replace the history file list for this day
        idt = idt0 + timedelta(days=nd)
        idt_str = datetime.strftime(idt,'%Y.%m.%d')
        print(' - working on ' + idt_str)
        sys.stdout.flush()
        fn_list = tfun.get_fn_list(idt, Ldir)
        
        # write the grid file (once per experiment) for plotting
        if write_grid == True:
            g_infile = fn_list[0]
            g_outfile = outdir / 'grid.nc'
            tfnc.write_grid(g_infile, g_outfile)
            write_grid = False

        # DO THE TRACKING
        if nd == 0: # first day
            # set IC
            plon0 = plon00.copy()
            plat0 = plat00.copy()
            pcs0 = pcs00.copy()
            # do the tracking
            if TR['start_hour'] > 0:
                fn_list = fn_list[TR['start_hour']:]
            P = tfun.get_tracks(fn_list, plon0, plat0, pcs0, TR, trim_loc=True)
            # save the results to NetCDF
            tfnc.start_outfile(out_fn, P)
        else: # subsequent days
            # set IC
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            # do the tracking
            P = tfun.get_tracks(fn_list, plon0, plat0, pcs0, TR)
            tfnc.append_to_outfile(out_fn, P)
        
    print(' - Took %0.1f sec for %s day(s)' %
            (time() - tt0, str(TR['days_to_track'])))
    print(50*'=')
print(50*'*' + '\nWrote to ' + str(outdir))

