"""
This creates and populates directories for ROMS runs on mox or similar.  It is
designed to work with the "BLANK" version of the .in file,
replacing things like $whatever$ with meaningful values.

To test from ipython on mac:

run make_dot_in -g cas6 -t v3 -x lo8b -r backfill -s continuation -d 2019.07.04 -bu 0 -np 196

"""

from pathlib import Path
import sys
from datetime import datetime
import pandas as pd

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import dot_in_argfun as dfun
import Lfun

Ldir = dfun.intro() # this handles all the argument passing

from datetime import datetime, timedelta
fdt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
fdt_yesterday = fdt - timedelta(1)

print('- dot_in.py creating files for LiveOcean for ' + Ldir['date_string'])

gtag = Ldir['gtag']
gtagex = Ldir['gtagex']
EX_NAME = Ldir['ex_name'].upper()

#### USER DEFINED VALUES ####

# which ROMS code to use
roms_name = 'LO_ROMS'

# account for differences when using biology
do_bio = True

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
    
ndtfast = 20
    
restart_nrrec = '-1' # '-1' for a non-crash restart file, otherwise '1' or '2'
his_interval = 3600 # seconds to define and write to history files
rst_interval = 10 # days between writing to the restart file (e.g. 5)

# which forcings to look for
atm_dir = 'atm1' # which atm forcing files to use
ocn_dir = 'ocn4' # which ocn forcing files to use
riv_dir = 'riv2' # which riv forcing files to use
tide_dir = 'tide2' # which tide forcing files to use

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

if dtsec == int(dtsec):
    dt = str(dtsec) + '.0d0' # a string version of dtsec, for the .in file
else:
    dt = str(dtsec) + 'd0' # a string version of dtsec, for the .in file
ninfo = int(his_interval/dtsec) # how often to write info to the log file (# of time steps)
nhis = int(his_interval/dtsec) # how often to write to the history files
ndefhis = int(nhis) # how often to create new history files
nrst = int(rst_interval*86400/dtsec)
ntimes = int(days_to_run*86400/dtsec)

# file location stuff
date_string = Ldir['date_string']
date_string_yesterday = fdt_yesterday.strftime(Lfun.ds_fmt)
dstart = str(int(Lfun.datetime_to_modtime(fdt) / 86400.))
f_string = 'f' + date_string
f_string_yesterday = 'f'+ date_string_yesterday

# Paths to forcing various file locations (all Path objects made into strings)
lo_dir = str(Ldir['LO'])
loo_dir = str(Ldir['LOo'])
grid_dir = str(Ldir['grid'])
Force_dir = Ldir['LOo'] / gtag / ('f' + Ldir['date_string'])
force_dir = str(Force_dir)
roms_dir = str(Ldir['roms'])

# hardwired because we don't have netCDF4
nrows = 1302 - 2
ncols = 663 - 2

# determine number of layers
S_fn = Ldir['grid'] / 'S_COORDINATE_INFO.csv'
S_info_dict = pd.read_csv(S_fn, index_col='ITEMS').to_dict()['VALUES']
nlayers = str(S_info_dict['N'])

# the .in file
dot_in_name = 'liveocean.in' # name of the .in file
out_dir = Ldir['roms'] / 'output' / Ldir['gtagex'] / ('f' + Ldir['date_string'])
out_dir_yesterday = Ldir['roms'] / 'output' / Ldir['gtagex'] / ('f' + date_string_yesterday)
Lfun.make_dir(out_dir, clean=True) # make sure it exists and is empty


if Ldir['start_type'] == 'continuation':
    nrrec = '0' # '-1' for a hot restart
    #ininame = 'ocean_rst.nc' # for a hot perfect restart
    ininame = 'ocean_his_0025.nc' # for a hot restart
    ini_fullname = str(out_dir_yesterday / ininame)
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini.nc' # could be an ini or history file
    ini_fullname = str(Force_dir / ocn_dir / ininame)

# END DERIVED VALUES

## create .in ##########################

f = open('BLANK.in','r')
f2 = open(out_dir / dot_in_name,'w')
in_varlist = ['base_dir','ntilei','ntilej','ntimes','dt','nrrec','ninfo',
    'nhis','dstart','ndefhis','nrst','force_dir','grid_dir','roms_dir',
    'atm_dir','ocn_dir','riv_dir','tide_dir',
    'ini_fullname','out_dir','EX_NAME','roms_name',
    'nrows','ncols', 'nlayers', 'ndtfast']
for line in f:
    for var in in_varlist:
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(eval(var)))
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
in_varlist = ['force_dir','riv_dir']
for line in f:
    for var in in_varlist:
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(eval(var)))
            line = line2
        else:
            line2 = line
    f3.write(line2)
f.close()
f3.close()
