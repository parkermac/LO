"""
This runs any of the forcing or post-processing jobs.

It can be run for a single forecast or over a range of past days.

- Test on mac from ipython:
run driver_forcing -g cas6 -t v3 -r backfill -s continuation -0 2019.07.04 -test True -f [FRC]
where [FRC] = ztest0, tide0, etc.

- Test on mac from command line:
python ./driver_forcing.py -g cas6 -t v3 -r backfill -s continuation -0 2019.07.04 -test True -f [FRC]

"""

import sys
import argparse
from datetime import datetime, timedelta
import subprocess

from lo_tools import Lfun, zfun

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas2
parser.add_argument('-t', '--tag', type=str)        # e.g. v3
parser.add_argument('-f', '--frc', type=str)        # e.g. tide
parser.add_argument('-r', '--run_type', type=str)   # forecast or backfill
parser.add_argument('-s', '--start_type', type=str) # new or continuation
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'tag', 'frc', 'run_type', 'start_type', 'ds0']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_forcing.py: ' + a)
        sys.exit()

if args.testing:
    from importlib import reload
    reload(Lfun)

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag)

# set time range to process
if args.run_type == 'forecast':
    ds0 = datetime.now().strftime(Lfun.ds_fmt)
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = dt0 + timedelta(days=Ldir['forecast_days'])
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
print((' Running %s for %s to %s ' % (args.run_type, ds0, ds1)).center(60,'-'))

# loop over all days
dt = dt0
while dt <= dt1:
    
    # make clean output directories
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + dt.strftime(Lfun.ds_fmt)) / args.frc
    print(('Creating %s' % (str(out_dir))).center(60,'-'))
    Lfun.make_dir(out_dir, clean=True)
    Lfun.make_dir(out_dir / 'Data')
    Lfun.make_dir(out_dir / 'Info')
    
    # run the code
    f_fn = Ldir['LO'] / 'forcing' / args.frc / 'make_forcing_main.py'
    cmd_list = ['python3', str(f_fn),
                '-g', args.gridname, '-t', args.tag, '-f', args.frc,
                '-r', args.run_type, '-s', args.start_type,
                '-d', dt.strftime(Lfun.ds_fmt), '-test', str(args.testing)]
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    with open(out_dir / 'Info' / 'screen_output.txt', 'w') as fout:
        fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(out_dir / 'Info' / 'subprocess_error.txt', 'w') as ffout:
            ffout.write(stderr.decode())
            
    if args.testing:
        print('\n' + ' sdtout '.center(60,'-'))
        print(stdout.decode())
        print('\n' + ' stderr '.center(60,'-'))
        print(stderr.decode())
        
    dt += timedelta(days=1)
    print('')
    



