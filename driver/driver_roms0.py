"""
This runs ROMS for one or more days, allowing for either a forecast or backfill.

Designed to be run from klone, and depends on other drivers having been run first on apogee.

The -np and -N flags specify the total number of cores, and the
cores per node.  For the current klone environment, acceptable choices are
-np any multiple of 40 up to 400, and -N 40
NOTE: np has to be an even multiple of N, and N has to be <= the number
of nodes of that size that I own.

To test on mac in ipython:
run driver_roms_mox -g cas6 -t v3t075 -x lo8 -r backfill -s continuation -0 2019.07.04 -np 200 -N 40 -v True --get_forcing False --run_roms False --move_his False

When running by hand on klone or mox it may help to use < /dev/null > test.log & at the end of the command.
The /dev/null input avoids occasionally having the job "stopped".

DEVELOPMENT NOTES: see the "various flags to facilitate testing" part of the arguments for other testing flags.
The -v (verbose) flag gives really useful screen output, so I have set its default to True.

2021.12.03 This version is copied from driver_roms.py and only differs in that the batch file
it creates is designed to work with our new LO_roms_source and LO_roms_user.

"""

import sys
import shutil
import argparse
from datetime import datetime, timedelta
from pathlib import Path
import subprocess
from time import time, sleep

# add the path by hand so that it will run on klone (outside of loenv)
pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas2
parser.add_argument('-t', '--tag', type=str)        # e.g. v3
parser.add_argument('-ta', '--tag_alt', type=str, default='') # used to make gtag in "remote_dir"
parser.add_argument('-x', '--ex_name', type=str)    # e.g. lo8b
parser.add_argument('-r', '--run_type', type=str)   # forecast or backfill
parser.add_argument('-s', '--start_type', type=str) # new or continuation
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-np', '--np_num', type=int) # e.g. 200, number of cores
parser.add_argument('-N', '--cores_per_node', type=int) # 40 on klone
# various flags to facilitate testing
parser.add_argument('-v', '--verbose', default=True, type=Lfun.boolean_string)
parser.add_argument('--get_forcing', default=True, type=Lfun.boolean_string)
parser.add_argument('--short_roms', default=False, type=Lfun.boolean_string)
parser.add_argument('--run_dot_in', default=True, type=Lfun.boolean_string)
parser.add_argument('--run_roms', default=True, type=Lfun.boolean_string)
parser.add_argument('--move_his', default=True, type=Lfun.boolean_string)
args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'np_num', 'cores_per_node']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_roms_klone.py: ' + a)
        sys.exit()
        
# set tag_alt to tag if it is not provided
if len(args.tag_alt) == 0:
    argsd['tag_alt'] = argsd['tag']

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)
Ldir['gtag_alt'] = Ldir['gridname'] + '_' + argsd['tag_alt']

# Assign some variables related to the remote machine.  These are user-specific
# and are specified in lo_tools/get_lo_info.py.
remote_user = Ldir['remote_user']
remote_machine = Ldir['remote_machine']
remote_dir0 = Ldir['remote_dir0']
local_user = Ldir['local_user']

# set time range to process
if args.run_type == 'forecast':
    ds0 = datetime.now().strftime(Lfun.ds_fmt)
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = dt0
    ds1 = dt1.strftime(Lfun.ds_fmt)
elif args.run_type == 'backfill': # you have to provide at least ds0 for backfill
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
print('Running ROMS %s %s-%s ' % (args.run_type, ds0, ds1))
sys.stdout.flush()

def messages(stdout, stderr, mtitle, test_flag):
    # utility function for displaying subprocess info
    if test_flag:
        print((' ' + mtitle + ' ').center(60,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
        sys.stdout.flush()

# Loop over days
dt = dt0
while dt <= dt1:
    f_string = 'f' + dt.strftime(Lfun.ds_fmt)
    print('>> ' + f_string + ' <<')
    print(' > started at %s' % (datetime.now().strftime('%Y.%m.%d %H:%M:%S')))
    sys.stdout.flush()
    
    # Set various paths.
    force_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / f_string
    dot_in_dir = Ldir['LO'] / 'dot_in' / Ldir['gtagex']
    dot_in_shared_dir = Ldir['LO'] / 'dot_in' / 'shared'
    roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    log_file = roms_out_dir / 'log.txt'
    #roms_ex_dir = Ldir['roms_code'] / 'makefiles' / Ldir['ex_name']
    roms_ex_dir = Ldir['parent'] / 'LO_roms_user' / Ldir['ex_name']
    if args.verbose:
        print(' - force_dir:    ' + str(force_dir))
        print(' - dot_in_dir:   ' + str(dot_in_dir))
        print(' - dot_in_shared_dir: ' + str(dot_in_shared_dir))
        print(' - roms_out_dir: ' + str(roms_out_dir))
        print(' - log_file:     ' + str(log_file))
        print(' - roms_ex_dir:  ' + str(roms_ex_dir))
        sys.stdout.flush()
    
    # Get the list of which forcing folders to copy
    force_dict = dict()
    with open(dot_in_dir / 'forcing_list.csv', 'r') as f:
        for line in f:
            which_force, force_choice = line.strip().split(',')
            force_dict[which_force] = force_choice
    
    if args.get_forcing:
        tt0 = time()
        # Name the place where the forcing files will be copied from
        remote_dir = remote_user + '@' + remote_machine + ':' + remote_dir0
        Lfun.make_dir(force_dir, clean=True)
        # Copy the forcing files, one folder at a time.
        for force in force_dict.keys():
            if force == 'open':
                pass
            else:
                force_choice = force_dict[force]
                cmd_list = ['scp','-r',
                    remote_dir + '/LO_output/forcing/' + Ldir['gtag_alt'] + '/' + f_string + '/' + force_choice,
                    str(force_dir)]
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                messages(stdout, stderr, 'Copy forcing ' + force_choice, args.verbose)
        print(' - time to get forcing = %d sec' % (time()-tt0))
        sys.stdout.flush()
    else:
        print(' ** skipped getting forcing')
        
    # Loop over blow ups.
    blow_ups = 0
    if args.short_roms:
        blow_ups_max = 0
    else:
        blow_ups_max = 5
    roms_worked = False
    while blow_ups <= blow_ups_max:
        print((' - Blow-ups = ' + str(blow_ups)))
        sys.stdout.flush()

        if args.run_dot_in:
            # Make the dot_in file.  NOTE: roms_out_dir is made clean by make_dot_in.py
            cmd_list = ['python3', str(dot_in_dir / 'make_dot_in.py'),
                        '-g', args.gridname, '-t', args.tag, '-x', args.ex_name,
                        '-r', args.run_type, '-s', args.start_type,
                        '-d', dt.strftime(Lfun.ds_fmt),
                        '-bu', str(blow_ups), '-np', str(args.np_num),
                        '-short_roms', str(args.short_roms)]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Make dot in', args.verbose)
            
            # Create batch script
            if 'klone' in Ldir['lo_env']:
                batch_name = 'make_klone_batch0.py'
            elif 'mox' in Ldir['lo_env']:
                batch_name = 'make_mox_batch.py'
            cmd_list = ['python3', str(dot_in_shared_dir / batch_name),
                '-xd', str(roms_ex_dir),
                '-rod', str(roms_out_dir),
                '-np', str(args.np_num),
                '-N', str(args.cores_per_node),
                '-x', Ldir['ex_name']]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Create batch script', args.verbose)
            
        else:
            print(' ** skipped making dot_in and batch script')
            args.run_roms = False # never run ROMS if we skipped making the dot_in
        
        if args.run_roms:
            tt0 = time()
            # Run ROMS using the batch script.
            if 'klone' in Ldir['lo_env']:
                cmd_list = ['sbatch', '-p', 'compute', '-A', 'macc','--wait',
                    str(roms_out_dir / 'klone_batch0.sh')]
            elif 'mox' in Ldir['lo_env']:
                cmd_list = ['sbatch', '-p', 'macc', '-A', 'macc','--wait',
                    str(roms_out_dir / 'mox_batch.sh')]
            # The --wait flag will cause the subprocess to not return until the job has terminated.
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Run ROMS', args.verbose)
            print(' - time to run ROMS = %d sec' % (time()-tt0))
            sys.stdout.flush()
    
            # A bit of checking to make sure that the log file exists...
            lcount = 0
            while not log_file.is_file():
                sleep(10)
                print(' - lcount = %d' % (lcount))
                sys.stdout.flush()
                lcount += 1
            # ...and that it is done being written to.
            llcount = 0
            log_done = False
            while log_done == False:
                sleep(3)
                if 'mox' in Ldir['lo_env']:
                    cmd_list = ['/usr/sbin/lsof', '-u', local_user,'|','grep',str(log_file)]
                else:
                    cmd_list = ['lsof', '-u', local_user,'|','grep',str(log_file)]
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                print(' - llcount = %d' % (lcount))
                sys.stdout.flush()
                lcount += 1
                if str(log_file) not in stdout.decode():
                    print(' - log done and closed')
                    sys.stdout.flush()
                    log_done = True
            # Look in the log file to see what happened, and decide what to do.
            roms_worked = False
            with open(log_file, 'r') as ff:
                found_line = False
                for line in ff:
                    if ('Blowing-up' in line) or ('BLOWUP' in line):
                        print(' - Run blew up, blow ups = ' + str(blow_ups))
                        found_line = True
                        roms_worked = False
                        break
                    elif 'ERROR' in line:
                        print(' - Run had an error. Check the log file.')
                        found_line = True
                        roms_worked = False
                        sys.exit()
                    elif 'ROMS/TOMS: DONE' in line:
                        found_line = True
                        print(' - ROMS SUCCESS')
                        roms_worked = True
                        break
            if not found_line:
                print(' - Problem finding line in log file.')
                sys.exit()
            sys.stdout.flush()
            if roms_worked:
                break # escape from blow_ups loop
            else:
                blow_ups += 1
        else:
            print(' ** skipped running ROMS')
            roms_worked = True
            break # escape from blow_ups loop

    if roms_worked:
        if args.move_his:
            tt0 = time()
            # Copy history files to the remote machine and clean up
            # (i) make sure the output directory exists
            cmd_list = ['ssh', remote_user + '@' + remote_machine,
                'mkdir -p ' + remote_dir0 + '/LO_roms/' + Ldir['gtagex']]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Make output directory on ' + remote_machine, args.verbose)
            # (ii) move the contents of roms_out_dir
            cmd_list = ['scp','-r',str(roms_out_dir),
                remote_user + '@' + remote_machine + ':' + remote_dir0 + '/LO_roms/' + Ldir['gtagex']]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Copy ROMS output to ' + remote_machine, args.verbose)
            # (iii) delete roms_out_dir and forcing files from the day before yesterday
            dt_prev = dt - timedelta(days=2)
            f_string_prev = 'f' + dt_prev.strftime(Lfun.ds_fmt)
            roms_out_dir_prev = Ldir['roms_out'] / Ldir['gtagex'] / f_string_prev
            force_dir_prev = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / f_string_prev
            shutil.rmtree(str(roms_out_dir_prev), ignore_errors=True)
            shutil.rmtree(str(force_dir_prev), ignore_errors=True)
            print(' - time to move history files and clean up = %d sec' % (time()-tt0))
            sys.stdout.flush()
        else:
            print(' ** skipped moving history files')
        dt += timedelta(days=1)
    else:
        print(' - ROMS FAIL for ' + dt.strftime(Lfun.ds_fmt))
        sys.exit()
        
    print(' > finished at %s' % (datetime.now().strftime('%Y.%m.%d %H:%M:%S')))
    sys.stdout.flush()
    



