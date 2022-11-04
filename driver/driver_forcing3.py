"""
This runs any of the forcing jobs.

It can be run for a single forecast or over a range of past days.

- Test on mac from ipython:
run driver_forcing3 -g cas6 -0 2019.07.04 -1 2019.07.05 -test True -f [FRC]
where [FRC] = riv00, etc.

- Test on mac from command line:
python driver_forcing3.py -g cas6 -0 2019.07.04 -test True -f [FRC]

As of 2022.06.08 the default start_type is perfect, which is the new default.
Not relevant in most cases. A start_type of "new" may be the only important use.

NEW: this "3" version uses the new organizational system where the forcing goes to a [gridname]
folder instead of [gtag].

"""

import sys
import argparse
from datetime import datetime, timedelta
import subprocess

from lo_tools import Lfun, zfun

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas2
parser.add_argument('-f', '--frc', type=str)        # e.g. tide
parser.add_argument('-r', '--run_type', type=str, default='backfill')   # forecast or backfill
parser.add_argument('-s', '--start_type', type=str, default='perfect') # new, continuation, or perfect
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'frc']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_forcing3.py: ' + a)
        sys.exit()

if args.testing:
    from importlib import reload
    reload(Lfun)

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname)

# set time range to process
if args.run_type == 'forecast':
    ds0 = datetime.now().strftime(Lfun.ds_fmt)
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = dt0
    ds1 = dt1.strftime(Lfun.ds_fmt)
elif args.run_type == 'backfill':
    ds0 = args.ds0
    if len(args.ds1) == 0:
        ds1 = ds0
    else:
        ds1 = args.ds1
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
else:
    print('Error: Unknown run_type')
    sys.exit()
    
if args.run_type == 'backfill':
    print(('Forcing: %s, %s, %s-%s ' % (args.frc, args.run_type, ds0, ds1)).center(60,'-'))
else:
    pass # don't clutter up the screen output for forecasts

# loop over all days
dt = dt0
while dt <= dt1:
    
    # fix start_type bug
    if (dt != dt0) and (args.start_type == 'new'):
        args.start_type = 'continuation'
    
    # make clean output directories
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + dt.strftime(Lfun.ds_fmt)) / args.frc
    Lfun.make_dir(out_dir, clean=True)
    Lfun.make_dir(out_dir / 'Data')
    Lfun.make_dir(out_dir / 'Info')
    
    # run the code
    f_fn = Ldir['LO'] / 'forcing' / args.frc / 'make_forcing_main.py'
    # If the user has this forcing in LO_user, then default to that:
    f_fn_alt = Ldir['LOu'] / 'forcing' / args.frc / 'make_forcing_main.py'
    if f_fn_alt.is_file():
        f_fn = f_fn_alt
        
    cmd_list = ['python3', str(f_fn),
                '-g', args.gridname, '-f', args.frc,
                '-r', args.run_type, '-s', args.start_type,
                '-d', dt.strftime(Lfun.ds_fmt), '-test', str(args.testing)]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    with open(out_dir / 'Info' / 'screen_output.txt', 'w') as fout:
        fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(out_dir / 'Info' / 'subprocess_error.txt', 'w') as ffout:
            ffout.write(stderr.decode())
    
    # this is intended to end up in the log that the cron job makes
    res_fn = out_dir / 'Info' / 'results.txt'
    if res_fn.is_file():
        with open(res_fn, 'r') as fout:
            for line in fout:
                print(line.replace('\n',''))
    else:
        print('ERROR: missing results.txt file')
            
    dt += timedelta(days=1)
    print('')
    



