"""
This runs ROMS for one or more days, allowing for either a forecast or backfill.

NOTES:

This version is intended to be an improvement on our previous drivers
(driver_roms[3,4,5,5a].py) with clearer handling of input parameters,
klone resource use, and error trapping.

It uses subprocess to watch for the end of a run, instead of pinging squeue.

It relies on LO/driver/batch/klone00_batch_BLANK.sh.

For run_type = forecast it assumes we are using the new (2025.08.26) scheme
in which foring always goes in individual day folders.

It sends the done_tags to their own folder LO/driver/done_tags

For testing/debugging these flags can be very useful:
-v True (verbose screen output)
--get_forcing False (don't get the forcing files, e.g. if they are already on klone)
--short_roms True (do a very short run, just a few timesteps)
--move_his False (don't move the results to apogee)

Testing on mac
run driver_roms00 -g cas7 -t t1 -x x11b -r backfill -0 2025.07.04 -s continuation -np 192 -cpu cpu-g2 -grp macc -vip True -v True --get_forcing False --run_roms False --move_his False

"""

import sys, os
import shutil
import argparse
from datetime import datetime, timedelta
from pathlib import Path
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time, sleep
import random
import string
from math import ceil

# Add the path to lo_tools by hand so that it we can import Lfun on klone
# without loenv. In general we write code to run on klone using only the
# default python3 installation.
pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

# >>> START Command Line Arguments <<<

# Arguments without defaults are required.

parser = argparse.ArgumentParser()
# Typically when you use these at the command line you can use the short version,
# like "-g" but the long version works as well "--gridname". In the code the longname
# is what is used.

# Basic info to specify which model configuration to run
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-t', '--tag', type=str)        # e.g. t1
parser.add_argument('-x', '--ex_name', type=str)    # e.g. x11b

# Set the run_type
parser.add_argument('-r', '--run_type', type=str, default='backfill')
# Choices: forecast or backfill

parser.add_argument('-s', '--start_type', type=str, default='continuation')
# Choices
# - new: only run 1 day, starting from ocean_ini.nc
# - perfect: start from ocean_rst.nc of previous day
# - continuation: start from ocean_his_[last one].nc of previous day
# - newperfect: first day = new, then perfect thereafter
# - newcontinuation: first day = new, then continuation thereafter

# If the run type is backfill you need to set the time range
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) #  this is set to ds0 if omitted

# The next two/three flags are needed for compute resource allocation.
#
# Set how many cpu's (cores) to use.
parser.add_argument('-np', '--np_num', type=int) # e.g. 160, number of cores
#
# Choose which type of computer resource
parser.add_argument('-cpu','--cpu_choice', type=str) # used in the sbatch script
# Choices: cpu-g2, compute, or ckpt-g2
#
# Choose which group (not needed when using --cpu_choice ckpt_g2)
parser.add_argument('-grp','--group_choice', type=str) # used in the sbatch script
# Choices: macc, coenv
# NOTE: Check with Parker before using macc or coenv.

# Optional flag used only for forecasts, so that they leave a clue if they are done for a given day.
# This allows a backup job in the crontab to exist but only run if needed. By using different
# done_tags we can have several different forecasts running, with backups, at the same time.
parser.add_argument('-done','--done_tag', type=str, default='00')

# Optional flags to facilitate testing.
parser.add_argument('-v', '--verbose', default=False, type=Lfun.boolean_string)
parser.add_argument('--get_forcing', default=True, type=Lfun.boolean_string)
parser.add_argument('--short_roms', default=False, type=Lfun.boolean_string)
parser.add_argument('--run_dot_in', default=True, type=Lfun.boolean_string)
parser.add_argument('--run_roms', default=True, type=Lfun.boolean_string)
parser.add_argument('--move_his', default=True, type=Lfun.boolean_string)

# Specialized flag, only for top-priority jobs like the daily forecasts 
parser.add_argument('-vip','--exclusive', default=False, type=Lfun.boolean_string)

# >>> END Command Line Arguments <<<

args = parser.parse_args()

# check for required arguments and other input errors
argsd = args.__dict__
for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'np_num', 'cpu_choice']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_roms##.py: ' + a)
        sys.exit()
if (argsd['cpu_choice'] == 'cpu-g2') and (argsd['group_choice'] not in ['macc','coenv']):
    print('*** Need a --group_choice for this --cpu_choice.')
    sys.exit()
if (argsd['cpu_choice'] == 'compute') and (argsd['group_choice'] != 'macc'):
    print('*** Need --group_choice macc for --cpu_choice compute.')
    sys.exit()
if (argsd['run_type'] == 'backfill') and (argsd['ds0'] == None):
    print('*** Need at least -0 YYYY.MM.DD for -r backfill')
    sys.exit()

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)

# Assign some variables related to the remote machine.  These are user-specific
# and are specified in get_lo_info.py.
remote_user = Ldir['remote_user']
remote_machine = Ldir['remote_machine']
remote_dir0 = Ldir['remote_dir0'] # NOTE: this is a string, not a Path object
local_user = Ldir['local_user']

# set time range to process
if args.run_type == 'forecast':
    ds0 = datetime.now().strftime(Lfun.ds_fmt)
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    # NOTE: to run as separate days we need Ldir['forecast_days'] to be an integer
    dt1 = dt0 + timedelta(days = (float(Ldir['forecast_days']) - 1))
    ds1 = dt1.strftime(Lfun.ds_fmt)
    # override for short_roms
    if args.short_roms == True:
        ds1 = ds0
        dt1 = dt0
    # Check to see if the forecast has already been run, and if so, exit!
    done_dir = Ldir['LO'] / 'driver' / 'done_tags'
    Lfun.make_dir(done_dir)
    done_fn = done_dir / ('forecast' + args.done_tag + '_done_' + ds0 + '.txt')
    if done_fn.is_file():
        print('Forecast has already run successfully - exiting')
        print(str(done_fn))
        sys.exit()
elif args.run_type == 'backfill': # you have to provide at least ds0 for backfill
    ds0 = args.ds0
    if args.ds1 == None:
        ds1 = ds0
    else:
        ds1 = args.ds1
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
else:
    print('Error: Unknown run_type')
    sys.exit()
print('Running ROMS %s %s-%s' % (args.run_type, ds0, ds1))
sys.stdout.flush()

Ncenter = 30
def messages(stdout, stderr, mtitle, verbose):
    # utility function for displaying subprocess info
    if verbose:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(Ncenter,'-'))
            print(stdout.decode())
    if len(stderr) > 0:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        # always print errors
        print(' stderr '.center(Ncenter,'-'))
        print(stderr.decode())
    sys.stdout.flush()

# RUNNING ROMS: Loop over days
dt = dt0
original_start_type = args.start_type
while dt <= dt1:
    
    # update start_type after the first day if needed
    
    if (dt == dt0) and (original_start_type in ['new','newperfect','newcontinuation']):
        start_type = 'new'
    else:
        start_type = original_start_type

    if (dt != dt0) and (original_start_type in ['newperfect','perfect']):
        start_type = 'perfect'

    if (dt != dt0) and (original_start_type in ['newcontinuation','continuation']):
        start_type = 'continuation'
        
    f_string = 'f' + dt.strftime(Lfun.ds_fmt)
    print('')
    print((' ' + f_string + ' ').center(Ncenter,'='))
    print(' > started at %s' % (datetime.now().strftime('%Y.%m.%d %H:%M:%S')))
    sys.stdout.flush()
    
    # Set various paths.
    force_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / f_string
    roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
    log_file = roms_out_dir / 'log.txt'
    
    # Look to see if there is a user instance of this dot_in, and if so, use it.
    user_dot_in_dir = Ldir['LOu'] / 'dot_in' / Ldir['gtagex']
    dot_in_dir = Ldir['LO'] / 'dot_in' / Ldir['gtagex']
    if user_dot_in_dir.is_dir():
        dot_in_dir = user_dot_in_dir
    
    # Decide where to look for ROMS executable.
    roms_ex_dir = Ldir['parent'] / 'LO_roms_user' / Ldir['ex_name']
    roms_ex_name = 'romsM'
    
    print(str(roms_out_dir)) # always print this
    if args.verbose:
        print(' - force_dir:    ' + str(force_dir))
        print(' - dot_in_dir:   ' + str(dot_in_dir))
        print(' - log_file:     ' + str(log_file))
        print(' - roms_ex_dir:  ' + str(roms_ex_dir))
        sys.stdout.flush()
    
    # Get the list of which forcing folders to copy
    force_dict = dict()
    with open(dot_in_dir / 'forcing_list.csv', 'r') as f:
        for line in f:
            which_force, force_choice = line.strip().split(',')
            force_dict[which_force] = force_choice
    
    # Get the forcing
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
                    remote_dir + '/LO_output/forcing/' + Ldir['gridname'] + '/' + f_string + '/' + force_choice,
                    str(force_dir)]
                proc = Po(cmd_list, stdout=Pi, stderr=Pi)
                stdout, stderr = proc.communicate()
                if len(stderr) > 0:
                    print('Error getting forcing %s' % (force))
                    sys.exit()
                    # got_forcing = False

                messages(stdout, stderr, 'Copy forcing ' + force_choice, args.verbose)
        print(' - time to get forcing = %d sec' % (time()-tt0))
    else:
        print(' ** skipped getting forcing')
        
    # Run ROMS. Loop over blow ups.
    blow_ups = 0
    if args.short_roms:
        blow_ups_max = 0
    else:
        blow_ups_max = 5
    roms_worked = False
    while blow_ups <= blow_ups_max:
        sys.stdout.flush()

        if args.run_dot_in:
            # Make the dot_in file.  NOTE: roms_out_dir is made clean by make_dot_in.py
            cmd_list = ['python3', str(dot_in_dir / 'make_dot_in.py'),
                        '-g', args.gridname, '-t', args.tag, '-x', args.ex_name,
                        '-r', 'backfill', '-s', start_type,
                        '-d', dt.strftime(Lfun.ds_fmt),
                        '-bu', str(blow_ups), '-np', str(args.np_num),
                        '-short_roms', str(args.short_roms)]
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Make dot in', args.verbose)
            
            # Create batch script
            f  = open(Ldir['LO'] / 'driver' / 'batch' / 'klone00_batch_BLANK.sh','r')
            f2 = open(roms_out_dir / 'klone_batch.sh','w')
            # generate a random jobname
            jobname = ''.join(random.choices(string.ascii_lowercase, k=5))

            # make choices about compute allocation
            np_num = args.np_num # integer
            if args.cpu_choice == 'compute': # compute nodes only have 40 cores
                ntasks_per_node = 40
                if np_num >= ntasks_per_node:
                    number_of_nodes = int(ceil(np_num / ntasks_per_node))
                else:
                    number_of_nodes = 1
                    ntasks_per_node = np_num
            if args.cpu_choice in ['cpu-g2', 'ckpt-g2']: # cpu-g2 nodes have 192 cores
                ntasks_per_node = 192
                if np_num >= ntasks_per_node:
                    number_of_nodes = int(ceil(np_num / ntasks_per_node))
                else:
                    number_of_nodes = 1
                    ntasks_per_node = np_num
            if args.exclusive:
                # only run this job on the node, and use all memory
                sbatch_exclusive_line = '#SBATCH --exclusive'
                sbatch_mem = '0'
            else:
                sbatch_exclusive_line = ''
                sbatch_mem = '128G'
                
            in_dict = {
                'jobname':jobname,
                'number_of_nodes':number_of_nodes,
                'ntasks_per_node':ntasks_per_node,
                'sbatch_mem':sbatch_mem,
                'sbatch_exclusive_line':sbatch_exclusive_line,
                'np_num':args.np_num,
                'roms_ex_dir':roms_ex_dir,
                'roms_ex_name':roms_ex_name,
                'roms_out_dir':roms_out_dir,
            }
            if args.verbose:
                print('Info for sbatch script:')
                for k in in_dict.keys():
                    print(' %s: %s' % (k, str(in_dict[k])))

            for line in f:
                for var in in_dict.keys():
                    if '$'+var+'$' in line: 
                        line2 = line.replace('$'+var+'$', str(in_dict[var]))
                        line = line2
                    else:
                        line2 = line
                f2.write(line2)
            f.close()
            f2.close()
        else:
            print(' ** skipped making dot_in and batch script')
            args.run_roms = False # never run ROMS if we skipped making the dot_in
        
        if args.run_roms:
            tt0 = time()
            # Run ROMS using the batch script.
            if args.group_choice != None:
                cmd_list = ['sbatch', '-p', args.cpu_choice, '-A', args.group_choice,
                    str(roms_out_dir / 'klone_batch.sh')]
            else:
                cmd_list = ['sbatch', '-p', args.cpu_choice,
                    str(roms_out_dir / 'klone_batch.sh')]
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                # Halt the program if there is a non-zero returncode from sbatch.
                print('EXITING due to sbatch error')
                print('sbatch returncode = %d' % (proc.returncode))
                sys.exit()
            messages(stdout, stderr, 'Run ROMS', args.verbose)
            print(' - time to run ROMS = %d sec' % (time()-tt0))
                    
            # Look in the log file to see what happened, and decide what to do.
            roms_worked = False
            with open(log_file, 'r') as ff:
                found_line = False
                for line in ff:
                    if ('Blowing-up' in line) or ('BLOWUP' in line) or ('Blows up' in line):
                        print(' - Run blew up, blow ups = ' + str(blow_ups))
                        found_line = True
                        roms_worked = False
                        
                        # save some info if the run blew up
                        print(' - blew up at %s' % (datetime.now().strftime('%Y.%m.%d %H:%M:%S')))
                        roms_bu_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / (f_string + '_blowup')
                        Lfun.make_dir(roms_bu_out_dir, clean=True)
                        try:
                            shutil.copy(roms_out_dir / 'log.txt', roms_bu_out_dir)
                            if args.verbose:
                                print(' - log.txt file saved to %s' % (str(roms_bu_out_dir)))
                        except FileNotFoundError:
                            print(' - log.txt file not found')
                        #
                        try:
                            his_list = roms_out_dir.glob('ocean_his_*.nc')
                            his_list = list(his_list)
                            his_list.sort()
                            if len(his_list) > 0:
                                shutil.copy(his_list[-1], roms_bu_out_dir)
                                if args.verbose:
                                    print(' - %s saved to %s' % (his_list[-1].name, str(roms_bu_out_dir)))
                            else:
                                print(' - no history files found')
                        except Exception as e:
                            print(' - problem saving blow-up history file:')
                            print(e)
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
                    elif 'ROMS: DONE' in line: # new version 2025.06.12
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
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Make output directory on ' + remote_machine, args.verbose)
            # (ii) move the contents of roms_out_dir
            cmd_list = ['scp','-r',str(roms_out_dir),
            remote_user + '@' + remote_machine + ':' + remote_dir0 + '/LO_roms/' + Ldir['gtagex']]
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Copy ROMS output to ' + remote_machine, args.verbose)
            # (iii) delete roms_out_dir and forcing files from several days in the past
            dt_prev = dt - timedelta(days=6)
            f_string_prev = 'f' + dt_prev.strftime(Lfun.ds_fmt)
            roms_out_dir_prev = Ldir['roms_out'] / Ldir['gtagex'] / f_string_prev
            roms_bu_out_dir_prev = Ldir['roms_out'] / Ldir['gtagex'] / (f_string_prev + '_blowup')
            force_dir_prev = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / f_string_prev
            shutil.rmtree(str(roms_out_dir_prev), ignore_errors=True)
            shutil.rmtree(str(roms_bu_out_dir_prev), ignore_errors=True)
            shutil.rmtree(str(force_dir_prev), ignore_errors=True)
            print(' - copy & clean up = %d sec' % (time()-tt0))
            sys.stdout.flush()
        else:
            print(' ** skipped moving history files')
            
        # Write a "done" file so that we know not to run the forecast again.
        if (dt == dt1) and (args.run_type == 'forecast'):
            with open(done_fn, 'w') as ffout:
                ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
            
        dt += timedelta(days=1)
    else:
        print(' - ROMS FAIL for ' + dt.strftime(Lfun.ds_fmt))
        sys.exit()
        
    print(' > finished at %s' % (datetime.now().strftime('%Y.%m.%d %H:%M:%S')))
    sys.stdout.flush()
    



