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

python3 driver_roms_mox.py -g cas6 -t v3 -x lo8b -r backfill -s continuation -0 2021.05.25 -np 196 -N 28 -test True > driver_log.txt &

or, after you have copied the forcing files once...

python3 driver_roms_mox.py -g cas6 -t v3 -x lo8b -r backfill -s continuation -0 2021.05.25 -np 196 -N 28 -test True -test2 True > driver_log.txt &

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
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
parser.add_argument('-test2', '--testing2', default=False, type=Lfun.boolean_string)
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

def messages(stdout, stderr, mtitle, test_flag):
    # utility function for displaying subprocess info
    if test_flag:
        print((' ' + mtitle + ' ').center(60,'='))
        print(' sdtout '.center(60,'-'))
        print(stdout.decode())
        print(' stderr '.center(60,'-'))
        print(stderr.decode())
        sys.stdout.flush()

# Loop over days
dt = dt0
while dt <= dt1:
    f_string = 'f' + dt.strftime(Lfun.ds_fmt)
    
    # Name the place where the forcing files will be copied from
    remote_dir='parker@boiler.ocean.washington.edu:/data1/parker'
    
    # Get the list of which forcing folders to copy
    dot_in_dir = Ldir['LO'] / 'dot_in' / Ldir['gtagex']
    force_dict = dict()
    with open(dot_in_dir / 'forcing_list.csv', 'r') as f:
        for line in f:
            which_force, force_choice = line.strip().split(',')
            force_dict[which_force] = force_choice
            
    # Make sure the directory exists where we are copying forcing files to.
    force_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / f_string
    
    if not args.testing2:
        Lfun.make_dir(force_dir, clean=True)
        # Copy the forcing files, one folder at a time.
        for force in force_dict.keys():
            force_choice = force_dict[force]
            cmd_list = ['scp','-r',
                remote_dir + '/LiveOcean_output/' + Ldir['gtag'] + '/' + f_string + '/' + force_choice,
                str(force_dir)]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Copy forcing ' + force_choice, args.testing)
    
    # Set some useful paths.
    roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    log_file = roms_out_dir / 'log.txt'
    roms_ex_dir = Ldir['roms_code'] / 'makefiles' / Ldir['ex_name']
    
    # Print helpful screen output.
    if args.testing:
        print(' - force_dir:    ' + str(force_dir))
        print(' - dot_in_dir:   ' + str(dot_in_dir))
        print(' - roms_out_dir: ' + str(roms_out_dir))
        print(' - roms_ex_dir:  ' + str(roms_ex_dir))
        print(' - log_file:     ' + str(log_file))
        sys.stdout.flush()
    
    # Loop over blow ups.
    blow_ups = 0
    if args.testing:
        blow_ups_max = 0
    else:
        blow_ups_max = 5
    roms_worked = False
    while blow_ups <= blow_ups_max:
        print((' Blow-ups = ' + str(blow_ups) + ' ').center(60,'.'))
        sys.stdout.flush()
    
        # Make the dot_in file.  NOTE: out_dir is made clean by make_dot_in.py
        f_fn = Ldir['LO'] / 'dot_in' / Ldir['gtagex'] / 'make_dot_in.py'
        cmd_list = ['python3', str(f_fn),
                    '-g', args.gridname, '-t', args.tag, '-x', args.ex_name,
                    '-r', args.run_type, '-s', args.start_type,
                    '-d', dt.strftime(Lfun.ds_fmt),
                    '-bu', str(blow_ups), '-np', str(args.np_num),
                    '-test', str(args.testing)]
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        messages(stdout, stderr, 'Make dot in', args.testing)
        
        # Create batch script
        cmd_list = ['python3', str(roms_ex_dir / 'make_back_batch_LO_version.py'),
            '-xp', str(roms_out_dir) +'/',
            '-np', str(args.np_num),
            '-N', str(args.cores_per_node),
            '-x', Ldir['ex_name']]
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        messages(stdout, stderr, 'Create batch script', args.testing)
            
        # Run ROMS using the batch script
        cmd_list = ['sbatch', '-p', 'macc', '-A', 'macc',
            str(roms_ex_dir / 'lo_back_batch_LO_version.sh')] # do not need , '&' when using subprocess?
        # ret1 = subprocess.call(cmd_list)
        # print('Return code = ' + str(ret1) + ' (0=success)')
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        messages(stdout, stderr, 'Run ROMS', args.testing)

        # Check the log file to see what happended, and decide the next step.
        # I think the subprocess for sbatch works as soon as the job is submitted, so we
        # need to keep looking for the log file until it exists.
        roms_worked = False
        sleep_sec = 10
        total_look_time = 5*3600
        max_look_count = int(total_look_time/sleep_sec)
        
        keep_looking = True
        look_count = 0
        flag = True
        while (look_count <= max_look_count) and keep_looking:
            print('-- Look count = ' + str(look_count))
            if log_file.is_file():
                if flag:
                    print('-- log file found')
                    flag = False
                with open(log_file, 'r') as ff:
                    for line in ff:
                        if ('Blowing-up' in line) or ('BLOWUP' in line):
                            print('Run blew up, blow ups = ' + str(blow_ups))
                            roms_worked = False
                            if args.testing:
                                print(line)
                            keep_looking = False
                            break
                        elif 'ERROR' in line:
                            print('Run had an error. Check the log file.')
                            roms_worked = False
                            if args.testing:
                                print(line)
                            sys.exit()
                        elif 'ROMS/TOMS: DONE' in line:
                            print('ROMS completed successfully.')
                            roms_worked = True
                            if args.testing:
                                print(line)
                            keep_looking = False
                            break
            if keep_looking:
                # it takes some time to write the log file, so even if we find it, we
                # may still have to keep checking for the text strings that allow us
                # to make a decision.
                time.sleep(sleep_sec)
                look_count += 1
            else:
                break # escape from while loop
        sys.stdout.flush()

        if roms_worked:
            break # escape from blow_ups loop
        else:
            blow_ups += 1

    if roms_worked:

        # TO DO: copy history files to boiler

        # TO DO: delete history files on mox for the day before yesterday

        dt += timedelta(days=1)
        
    else:
        print('ROMS did not work for ' + dt.strftime(Lfun.ds_fmt))
        sys.exit()
    sys.stdout.flush()
    



