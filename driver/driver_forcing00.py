"""
This runs any of the forcing jobs.

It can be run for a single forecast or over a range of past days.

Test on mac from ipython:
run driver_forcing00 -g [gridname] -0 2019.07.04 -1 2019.07.05 -test True -f [FRC]
where [FRC] = atm01, etc.

NOTE: This is just like driver_forcing3.py except it does not use clean=True
when creating the output directory for a day IF the run_type = forecast.
This change is because we don't want to lose previous forcing if a given
forecast is going to planB.

There is also a new optional input: -test_planB True to test planB.

This this driver is suitable for the new forcing organization scheme
(August 2025) where everything goes in single day folders, including forecasts.
Thus for a forecast, we now populate three separate day folders instead of one.

Because of this change you should ONLY use this driver for forcing which has
this new forecast organization scheme implemented.

Current examples of compliant forcing are: atm01, ocnG01, riv01, tide01, trapsF??
and ocnN.

This is meant to then be paired with driver_roms00.py which will expect forecasts
to be organized this way.

"""

import sys
import argparse
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun, zfun

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-f', '--frc', type=str)        # e.g. tide01
parser.add_argument('-r', '--run_type', type=str, default='backfill')   # forecast or backfill
parser.add_argument('-s', '--start_type', type=str, default='continuation') # new, continuation, or perfect
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
# new optional argument to test whatever planB is implemented for a given forcing
parser.add_argument('-test_planB', default=False, type=zfun.boolean_string)

# optional arguments used only for ocnN, to determine what to nest inside
parser.add_argument('-gtx', '--gtagex', type=str) # e.g. cas7_t1_x11b
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 0 = Ldir['roms_out1'], etc.
parser.add_argument('-do_bio', default=False, type=Lfun.boolean_string) # True to add bio vars to ocnN forcing

# optional argument to select different version of traps pre
parser.add_argument('-tP', '--trapsP', type=str, default='trapsP00') # LO/pre/trapsP## version

args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'frc']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_forcing00.py: ' + a)
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
    
    # make output directories
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
        ('f' + dt.strftime(Lfun.ds_fmt)) / args.frc
    if args.run_type == 'backfill':
        Lfun.make_dir(out_dir, clean=True)
    elif args.run_type == 'forecast':
        Lfun.make_dir(out_dir)
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
                '-tP', args.trapsP,
                '-d', dt.strftime(Lfun.ds_fmt), '-test', str(args.testing),
                '-test_planB', str(args.test_planB),
                '-gtx', args.gtagex, '-ro', str(args.roms_out_num), '-do_bio', str(args.do_bio)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
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
    sys.stdout.flush()
            
    dt += timedelta(days=1)
    print('')
    



