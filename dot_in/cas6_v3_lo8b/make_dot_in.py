"""
This creates and populates directories for ROMS runs on mox or similar.  It is
designed to work with the "BLANK" version of the .in file,
replacing things like $whatever$ with meaningful values.

To test from ipython on mac:
run make_dot_in -g cas6 -t v3 -x lo8b -r backfill -s continuation -d 2019.07.04 -bu 0 -np 196

If you are using an alternate tag then run like:
run make_dot_in -g cas6 -t v3 -t v3a -x lo8b -r backfill -s continuation -d 2019.07.04 -bu 0 -np 196
But note that the resulting gtagex should be the same as the folder this is in, and it
will throw a WARNING pointing this out if it is not.
"""

# NOTE: we limit the imports to modules that exist in python3 on mox
from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import dot_in_argfun as dfun
import Lfun

Ldir = dfun.intro() # this handles all the argument passing

# check that the folder name matches the arguments (cludgey)
ppth = Path(__file__).absolute().parent.name
gridname_check, tag_check, ex_name_check = (str(ppth)).split('_')
if Ldir['gridname'] != gridname_check:
    print('WARNING: gridname mismatch')
if Ldir['tag'] != tag_check:
    print('WARNING: tag mismatch')
if Ldir['ex_name'] != ex_name_check:
    print('WARNING: ex_name mismatch')

from datetime import datetime, timedelta
fdt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
fdt_yesterday = fdt - timedelta(1)

print('Creating dot_in for ' + Ldir['date_string'])
# initialize dict to hold values that we will substitute into the dot_in file.
D = dict()

D['EX_NAME'] = Ldir['ex_name'].upper()

#### USER DEFINED VALUES ####

# which ROMS code to use
D['roms_name'] = 'LO_ROMS'

multi_core = True # use more than one core

if Ldir['run_type'] == 'backfill':
    days_to_run = 1.0
elif Ldir['run_type'] == 'forecast':
    days_to_run = float(Ldir['forecast_days'])

# time step in seconds (should fit evenly into 3600 sec)
if Ldir['blow_ups'] == 0:
    dtsec = 40
elif Ldir['blow_ups'] == 1:
    dtsec = 30
elif Ldir['blow_ups'] == 2:
    dtsec = 25
elif Ldir['blow_ups'] == 3:
    dtsec = 20
elif Ldir['blow_ups'] == 4:
    dtsec = 15
elif Ldir['blow_ups'] == 5:
    dtsec = 10
else:
    print('Unsupported number of blow ups: %d' % (Ldir['blow_ups']))

D['ndtfast'] = 20
    
his_interval = 3600 # seconds to define and write to history files
rst_interval = 10 # days between writing to the restart file (e.g. 5)

# which forcings to look for
D['atm'] = 'atm1' # which atm forcing files to use
D['ocn'] = 'ocn4' # which ocn forcing files to use
D['riv'] = 'riv2' # which riv forcing files to use
D['tide'] = 'tide2' # which tide forcing files to use

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
    elif Ldir['np_num'] == 392:
        ntilei = '14' # number of tiles in I-direction
        ntilej = '28' # number of tiles in J-direction
    elif Ldir['np_num'] == 588:
        ntilei = '21' # number of tiles in I-direction
        ntilej = '28' # number of tiles in J-direction
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

D['ninfo'] = int(his_interval/dtsec) # how often to write info to the log file (# of time steps)
D['nhis'] = int(his_interval/dtsec) # how often to write to the history files
D['ndefhis'] = D['nhis'] # how often to create new history files
D['nrst'] = int(rst_interval*86400/dtsec)
D['ntimes'] = int(days_to_run*86400/dtsec)

# file location stuff
date_string_yesterday = fdt_yesterday.strftime(Lfun.ds_fmt)
D['dstart'] = int(Lfun.datetime_to_modtime(fdt) / 86400.)

# Paths to forcing various file locations
D['grid_dir'] = Ldir['grid']
force_dir = Ldir['LOo'] / Ldir['gtag'] / ('f' + Ldir['date_string'])
D['force_dir'] = force_dir
D['roms_dir'] = Ldir['roms']

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
out_dir = Ldir['roms'] / 'output' / Ldir['gtagex_alt'] / ('f' + Ldir['date_string'])
D['out_dir'] = out_dir
out_dir_yesterday = Ldir['roms'] / 'output' / Ldir['gtagex_alt'] / ('f' + date_string_yesterday)
Lfun.make_dir(out_dir, clean=True) # make sure it exists and is empty

if Ldir['start_type'] == 'continuation':
    nrrec = '0' # '-1' for a hot restart
    #ininame = 'ocean_rst.nc' # for a hot perfect restart
    ininame = 'ocean_his_0025.nc' # for a hot restart
    ini_fullname = out_dir_yesterday / ininame
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini.nc' # could be an ini or history file
    ini_fullname = force_dir / D['ocn'] / ininame
D['nrrec'] = nrrec
D['ini_fullname'] = ini_fullname

# END DERIVED VALUES

## create .in ##########################
f = open('BLANK.in','r')
f2 = open(out_dir / 'liveocean.in','w')
for line in f:
    for var in D.keys():
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(D[var]))
            line = line2
        else:
            line2 = line
    f2.write(line2)
f.close()
f2.close()

## npzd2o_Banas.in ###########

f = open('npzd2o_Banas_BLANK.in','r')
bio_dot_in_name = 'npzd2o_Banas.in'
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
