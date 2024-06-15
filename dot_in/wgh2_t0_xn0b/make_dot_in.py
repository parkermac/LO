"""
This creates and populates directories for ROMS runs on mox or similar.
Its main product is the .in file used for a ROMS run.

It is designed to work with the BLANK.in file, replacing things like
$whatever$ with meaningful values.

To test from ipython on mac:
run make_dot_in -g wgh0 -t t0 -x xn0 -r backfill -s continuation -d 2021.07.04 -bu 0 -np 64

If you call with -short_roms True it will create a .in that runs for a shorter time or
writes history files more frequently (exact behavior is in the code below).  This can
be really useful for debugging.

"""

# NOTE: we limit the imports to modules that exist in python3 on mox
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import dot_in_argfun as dfun
import Lfun

dot_in_dir = Path(__file__).absolute().parent

Ldir = dfun.intro() # this handles all the argument passing

# check that the folder name matches the arguments (cludgey)
ppth = Path(__file__).absolute().parent.name
gridname_check, tag_check, ex_name_check = (str(ppth)).split('_')
if Ldir['gridname'] != gridname_check:
    print('WARNING: gridname mismatch')
    sys.exit()
if Ldir['tag'] != tag_check:
    print('WARNING: tag mismatch')
    sys.exit()
if Ldir['ex_name'] != ex_name_check:
    print('WARNING: ex_name mismatch')
    sys.exit()

fdt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
fdt_yesterday = fdt - timedelta(days=1)

print(' --- making dot_in for ' + Ldir['date_string'])
# initialize dict to hold values that we will substitute into the dot_in file.
D = dict()

D['EX_NAME'] = Ldir['ex_name'].upper()

#### USER DEFINED VALUES ####

# which ROMS code to use
D['roms_name'] = 'LO_roms_source'

multi_core = True # use more than one core

if Ldir['run_type'] == 'backfill':
    days_to_run = 1.0
elif Ldir['run_type'] == 'forecast':
    days_to_run = float(Ldir['forecast_days'])

# time step in seconds (should fit evenly into 3600 sec)
if Ldir['blow_ups'] == 0:
    dtsec = 12
elif Ldir['blow_ups'] == 1:
    dtsec = 9
elif Ldir['blow_ups'] == 2:
    dtsec = 6
elif Ldir['blow_ups'] == 3:
    dtsec = 4
elif Ldir['blow_ups'] == 4:
    dtsec = 3
elif Ldir['blow_ups'] == 5:
    dtsec = 2
else:
    print('Unsupported number of blow ups: %d' % (Ldir['blow_ups']))

D['ndtfast'] = 20
    
his_interval = 3600 # seconds to define and write to history files
rst_interval = 1 # days between writing to the restart file (e.g. 5)

# Find which forcings to look for (search the csv file in this directory).
# We use the csv file because driver_roms_mox.py also uses it to copy forcing
# without extra stuff.
this_dir = Path(__file__).absolute().parent
with open(this_dir / 'forcing_list.csv', 'r') as f:
    for line in f:
        which_force, force_choice = line.strip().split(',')
        D[which_force] = force_choice
        # Note that the last of these will not be a force, but will
        # instead encode the open boundaries, e.g.: D['open'] = 'NSW'
        
# populate open boundary choices
for O in list('NSEW'):
    if O in list(D['open']):
        D[O+'fs'] = 'Cha'
        D[O+'2']  = 'Fla'
        D[O+'3']  = 'RadNud'
    else:
        D[O+'fs'] = 'Clo'
        D[O+'2']  = 'Clo'
        D[O+'3']  = 'Clo'

#### END USER DEFINED VALUES ####

# DERIVED VALUES
if multi_core:
    if Ldir['np_num'] == 64: # for new mox nodes 2*32=64 2019_02
        ntilei = '8' # number of tiles in I-direction
        ntilej = '8' # number of tiles in J-direction
    elif Ldir['np_num'] == 72:
        ntilei = '6' # number of tiles in I-direction
        ntilej = '12' # number of tiles in J-direction
    elif Ldir['np_num'] == 112:
        ntilei = '8' # number of tiles in I-direction
        ntilej = '14' # number of tiles in J-direction
    elif Ldir['np_num'] == 144:
        ntilei = '8' # number of tiles in I-direction
        ntilej = '18' # number of tiles in J-direction
    elif Ldir['np_num'] == 196:
        ntilei = '14' # number of tiles in I-direction
        ntilej = '14' # number of tiles in J-direction
    elif Ldir['np_num'] == 252:
        ntilei = '14' # number of tiles in I-direction
        ntilej = '18' # number of tiles in J-direction
    elif Ldir['np_num'] == 392:
        ntilei = '14' # number of tiles in I-direction
        ntilej = '28' # number of tiles in J-direction
    elif Ldir['np_num'] == 588:
        ntilei = '21' # number of tiles in I-direction
        ntilej = '28' # number of tiles in J-direction
    elif Ldir['np_num'] == 400: # klone
        ntilei = '20' # number of tiles in I-direction
        ntilej = '20' # number of tiles in J-direction
    elif Ldir['np_num'] == 200: # klone
        ntilei = '10' # number of tiles in I-direction
        ntilej = '20' # number of tiles in J-direction
    elif Ldir['np_num'] == 40: # klone
        ntilei = '5' # number of tiles in I-direction
        ntilej = '8' # number of tiles in J-direction
    elif Ldir['np_num'] == 4: # klone
        ntilei = '2' # number of tiles in I-direction
        ntilej = '2' # number of tiles in J-direction
    elif Ldir['np_num'] == 160: # klone g2 nodes (have 5*32 = 160)
        ntilei = '10' # number of tiles in I-direction
        ntilej = '16' # number of tiles in J-direction
    else:
        print('Unsupported number of processors: %d' % (Ldir['np_num']))
else:
    ntilei = '1'
    ntilej = '1'
D['ntilei'] = ntilei
D['ntilej'] = ntilej

# a string version of dtsec, for the .in file
if dtsec == int(dtsec):
    dt = str(dtsec) + '.0d0'
else:
    dt = str(dtsec) + 'd0'
D['dt'] = dt

if Ldir['short_roms']:
    print(' --- running short roms')
    his_interval = 10 * dtsec
    D['ntimes'] = int(his_interval/dtsec) # run for his_interval
else:
    D['ntimes'] = int(days_to_run*86400/dtsec)

D['ninfo'] = int(his_interval/dtsec) # how often to write info to the log file (# of time steps)
D['nhis'] = int(his_interval/dtsec) # how often to write to the history files
D['ndefhis'] = D['nhis'] # how often to create new history files
D['nrst'] = int(rst_interval*86400/dtsec)

# file location stuff
date_string_yesterday = fdt_yesterday.strftime(Lfun.ds_fmt)
D['dstart'] = int(Lfun.datetime_to_modtime(fdt) / 86400.)

# Paths to forcing various file locations
D['grid_dir'] = Ldir['grid']
force_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string'])
D['force_dir'] = force_dir
D['roms_varinfo_dir'] = Ldir['parent'] / 'LO_roms_source_alt' / 'varinfo'

# get horizontal coordinate info
with open(Ldir['grid'] / 'XY_COORDINATE_INFO.csv','r') as xyf:
    for line in xyf:
        ltup = line.split(',')
        D[ltup[0]] = int(ltup[1]) - 2

# get vertical coordinate info
with open(Ldir['grid'] / 'S_COORDINATE_INFO.csv','r') as sf:
    for line in sf:
        ltup = line.split(',')
        if ltup[0] != 'ITEMS':
            D[ltup[0]] = int(ltup[1])

# the output directory and the one from the day before
out_dir = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string'])
D['out_dir'] = out_dir
out_dir_yesterday = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + date_string_yesterday)
Lfun.make_dir(out_dir, clean=True) # make sure it exists and is empty

if Ldir['start_type'] == 'perfect':
    nrrec = '-1' # '-1' for a perfect restart
    ininame = 'ocean_rst.nc' # start from restart file
    ini_fullname = out_dir_yesterday / ininame
elif Ldir['start_type'] == 'continuation':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_his_0025.nc' # restart from a history file
    ini_fullname = out_dir_yesterday / ininame
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini.nc' # start from ini file
    ini_fullname = force_dir / D['ocn'] / ininame
D['nrrec'] = nrrec
D['ini_fullname'] = ini_fullname

# END DERIVED VALUES

## create liveocean.in ##########################
f = open(dot_in_dir / 'BLANK.in','r')
f2 = open(out_dir / 'liveocean.in','w')
for line in f:
    for var in D.keys():
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(D[var]))
            line = line2 # needed because we loop over all "var" for line
        else:
            line2 = line
    f2.write(line2)
f.close()
f2.close()

# create bio_Banas.in #######################
f = open(dot_in_dir / 'bio_Fennel_BLANK.in','r')
bio_dot_in_name = 'bio_Banas.in'
f3 = open(out_dir / bio_dot_in_name,'w')
for line in f:
    for var in D.keys():
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(D[var]))
            line = line2
        else:
            line2 = line
    f3.write(line2)
f.close()
f3.close()
