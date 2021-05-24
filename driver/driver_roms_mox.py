"""
This runs ROMS for one or more days, allowing for either a forecast or backfill.

Designed to be run from mox, and depends on other drivers having been run first on boiler

The -np and -N flags specify the total number of cores, and the
cores per node.  For the current mox environment, acceptable choices are
-np 64 -N 32
  or
-np 196 -N 28
i.e. np has to be an even multiple of N, and N has to be <= the number
of nodes of that size that I own.

To test on mac in ipython

run driver_roms_mox -g cas6 -t v3 -x lo8b -r backfill -s continuation -0 2019.07.04 -np 196 -N 28

to test on mox

python3 driver_roms_mox.py -g cas6 -t v3 -x lo8b -r backfill -s continuation -0 2021.05.24 -np 196 -N 28 -test True

"""

import sys
import argparse
from datetime import datetime, timedelta
from pathlib import Path
import subprocess
import time

pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

def boolean_string(s):
    # used by argparse
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas2
parser.add_argument('-t', '--tag', type=str)        # e.g. v3
parser.add_argument('-x', '--ex_name', type=str)    # e.g. lo8b
parser.add_argument('-f', '--frc', type=str)        # e.g. tide
parser.add_argument('-r', '--run_type', type=str)   # forecast or backfill
parser.add_argument('-s', '--start_type', type=str) # new or continuation
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-np', '--np_num', type=int) # e.g. 196, number of cores
parser.add_argument('-N', '--cores_per_node', type=int) # 28 or 32 on mox, number of cores per node
parser.add_argument('-test', '--testing', default=False, type=boolean_string)
args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'ds0', 'np_num', 'cores_per_node']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_roms_mox.py: ' + a)
        sys.exit()

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)

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
print((' Running ROMS %s for %s to %s ' % (args.run_type, ds0, ds1)).center(60,'*'))
sys.stdout.flush()

# Loop over days
dt = dt0
while dt <= dt1:
    
    # Copy the forcing files for this day from the computer that made them (e.g. boiler)
    remote_dir='parker@boiler.ocean.washington.edu:/data1/parker'
    f_string = 'f' + dt.strftime(Lfun.ds_fmt)
    # make sure the directory exists where we are copying files to
    Lfun.make_dir(Ldir['LOo'] / 'forcing' / Ldir['gtag'])
    
    if not args.testing:
        cmd_list = ['scp','-r',
            remote_dir + '/LiveOcean_output/' + Ldir['gtag'] + '/' + f_string,
            str(Ldir['LOo']) + '/forcing/' + Ldir['gtag'] + '/' + f_string]
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if args.testing:
            print(' Copy forcing '.center(60,'='))
            print('\n' + ' sdtout '.center(60,'-'))
            print(stdout.decode())
            print('\n' + ' stderr '.center(60,'-'))
            print(stderr.decode())
            sys.stdout.flush()
        
    
    # Set some useful paths
    roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    log_file = roms_out_dir / 'log.txt'
    roms_ex_dir = Ldir['roms_code'] / 'makefiles' / Ldir['ex_name']
    
    # Loop over blow ups
    blow_ups = 0
    if args.testing:
        blow_ups_max = 0
    else:
        blow_ups_max = 7
        
    roms_worked = False
    while blow_ups <= blow_ups_max:
    
        # Make the dot_in file.  NOTE: out_dir is made clean by make_dot_in.py
        f_fn = Ldir['LO'] / 'dot_in' / Ldir['gtagex'] / 'make_dot_in.py'
        cmd_list = ['python3', str(f_fn),
                    '-g', args.gridname, '-t', args.tag, '-x', args.ex_name,
                    '-r', args.run_type, '-s', args.start_type,
                    '-d', dt.strftime(Lfun.ds_fmt),
                    '-bu', str(blow_ups), '-np', str(args.np_num)]
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if args.testing:
            print(' Make dot in '.center(60,'='))
            print('\n' + ' sdtout '.center(60,'-'))
            print(stdout.decode())
            print('\n' + ' stderr '.center(60,'-'))
            print(stderr.decode())
            sys.stdout.flush()
        
        # Create batch script
        cmd_list = ['python3', str(roms_ex_dir / 'make_back_batch_LO_version.py'),
            '-xp', str(roms_out_dir) +'/',
            '-np', str(args.np_num),
            '-N', str(args.cores_per_node),
            '-x', Ldir['ex_name']]
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if args.testing:
            print(' Create batch script '.center(60,'='))
            print('\n' + ' sdtout '.center(60,'-'))
            print(stdout.decode())
            print('\n' + ' stderr '.center(60,'-'))
            print(stderr.decode())
            sys.stdout.flush()
            

    #     # Run ROMS using the batch script
    #     cmd_list = ['sbatch', '-p', 'macc', '-A', 'macc',
    #         str(roms_ex_dir / 'lo_back_batch_LO_version.sh'), '&']
    #     ret1 = subprocess.call(cmd_list)
    #     print('Return code = ' + str(ret1) + ' (0=success)')
    #
    #     # check log file to see if it worked
    #     roms_worked = False
    #     keep_looking = True
    #     while keep_looking:
    #         time.sleep(10)
    #         if log_file.is_file():
    #             with open(log_file, 'r') as ff:
    #                 for line in ff:
    #                     if ('Blowing-up' in line) or ('BLOWUP' in line):
    #                         roms_worked = False
    #                         keep_looking = False
    #                         break
    #                     elif 'ERROR' in line:
    #                         roms_worked = False
    #                         keep_looking = False
    #                         break
    #                     elif 'ROMS/TOMS: DONE' in line:
    #                         roms_worked = True
    #                         keep_looking = False
    #                         break
    #
        if args.testing:
            roms_worked = True
            
        if roms_worked: # test that run completed successfully
            break
        else:
            blow_ups += 1
    #
    if roms_worked:

        # TO DO: copy history files to boiler

        # TO DO: delete history files on mox for the day before yesterday

        dt += timedelta(days=1)
    # else:
    #     print('ROMS did not work for ' + dt.strftime(Lfun.ds_fmt))
    #     sys.exit()
    



