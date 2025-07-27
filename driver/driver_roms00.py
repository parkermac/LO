"""
This runs ROMS for one or more days, allowing for either a forecast or backfill.

This version is intended to be an improvement on our previous drivers
(driver_roms[3,4,5,5a].py) with clearer handling of input parameters,
klone resource use, and error trapping. In particular I am experimenting
with using subprocess to watch for the end of a run, instead of
pinging squeue.

It relies on LO/driver/batch/klone00_make_batch.py and klone00_batch_BLANK.sh.

For testing/debugging these flags can be very useful:
-v True (verbose screen output)
--get_forcing False (don't get the forcing files, e.g. if they are already on klone)
--short_roms True (do a very short run, just a few timesteps)
--move_his False (don't move the results to apogee)

"""

import sys, os
import shutil
import argparse
from datetime import datetime, timedelta
from pathlib import Path
import subprocess
from time import time, sleep
import random
import string

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
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted

# The next three flags are needed for compute resource allocation.

# Specify how many cpu's (cores) to use
parser.add_argument('-tpn', '--tasks_per_node', type=int) # e.g. 160, number of cores to use, up to 192

# Optional flag to use a different type of node/slice
parser.add_argument('--cpu_choice', type=str, default='cpu-g2') # used in the sbatch call
# Generally use cpu-gs unless you specifically want to use our older klone nodes.
# Choices: cpu-g2, compute, or ckpt-g2

# Optional flag to use a different group choice.
parser.add_argument('--group_choice', type=str, default='macc') # used in the sbatch call
# Choices: macc, coenv
# NOTE: Check with Parker before using macc or coenv. Not used when cpu_choice = ckpt-g2.

# Optional flag used only for forecasts, so that they leave a clue if they are done for a given day.
# This allows a backup job in the crontab to exist but only run if needed.
parser.add_argument('--done_tag', type=str, default='3') # used in done_fn

# various flags to facilitate testing
parser.add_argument('-v', '--verbose', default=False, type=Lfun.boolean_string)
parser.add_argument('--get_forcing', default=True, type=Lfun.boolean_string)
parser.add_argument('--short_roms', default=False, type=Lfun.boolean_string)
parser.add_argument('--run_dot_in', default=True, type=Lfun.boolean_string)
parser.add_argument('--run_roms', default=True, type=Lfun.boolean_string)
parser.add_argument('--move_his', default=True, type=Lfun.boolean_string)

# >>> END Command Line Arguments <<<

args = parser.parse_args()

# check for required arguments
argsd = args.__dict__
for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'tasks_per_node']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_roms3.py: ' + a)
        sys.exit()

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)

# Assign some variables related to the remote machine.  These are user-specific
# and are specified in get_lo_info.py.
remote_user = Ldir['remote_user']
remote_machine = Ldir['remote_machine']
remote_dir0 = Ldir['remote_dir0']
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
    done_fn = Ldir['LO'] / 'driver' / ('forecast' + args.done_tag + '_done_' + ds0 + '.txt')
    if done_fn.is_file():
        print('Forecast has already run successfully - exiting')
        print(str(done_fn))
        sys.exit()
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

# Loop over days
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
    f_string0 = 'f' + dt0.strftime(Lfun.ds_fmt) # used for forcing a forecast
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
    
    # Decide where to look for executable.
    # use updated ROMS
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
    
    if args.get_forcing:
        for fff in range(10):
            # We put this in a loop to allow it to try several times. This is prompted
            # by intermittent ssh_exchange_identification errors, particularly on mox.
            got_forcing = True
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
                    if (args.run_type == 'backfill') or (force_choice == 'ocnN'):
                        F_string = f_string
                    elif args.run_type == 'forecast':
                        F_string = f_string0
                    cmd_list = ['scp','-r',
                        remote_dir + '/LO_output/forcing/' + Ldir['gridname'] + '/' + F_string + '/' + force_choice,
                        str(force_dir)]
                    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = proc.communicate()
                    if len(stderr) > 0:
                        got_forcing = False
                    messages(stdout, stderr, 'Copy forcing ' + force_choice, args.verbose)
            print(' - time to get forcing = %d sec' % (time()-tt0))
            sys.stdout.flush()
            if got_forcing == True:
                break
            else:
                sleep(60)
        if got_forcing == False:
            print('Error getting forcing, fff = %d' % (fff))
            sys.exit()
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
        # print((' - Blow-ups = ' + str(blow_ups)))
        sys.stdout.flush()

        if args.run_dot_in:
            # Make the dot_in file.  NOTE: roms_out_dir is made clean by make_dot_in.py
            cmd_list = ['python3', str(dot_in_dir / 'make_dot_in.py'),
                        '-g', args.gridname, '-t', args.tag, '-x', args.ex_name,
                        '-r', 'backfill', '-s', start_type,
                        '-d', dt.strftime(Lfun.ds_fmt),
                        '-bu', str(blow_ups), '-np', str(args.tasks_per_node),
                        '-short_roms', str(args.short_roms)]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Make dot in', args.verbose)
            
            # Create batch script
            batch_name = 'klone5_make_batch.py'
            # generate a random jobname
            jobname = ''.join(random.choices(string.ascii_lowercase, k=5))
            cmd_list = ['python3', str(Ldir['LO'] / 'driver' / 'batch' / batch_name),
                '-xd', str(roms_ex_dir),
                '-rxn', roms_ex_name,
                '-rod', str(roms_out_dir),
                '-tpn', str(args.tasks_per_node),
                '-x', Ldir['ex_name'],
                '-j', jobname]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            messages(stdout, stderr, 'Create batch script', args.verbose)
        else:
            print(' ** skipped making dot_in and batch script')
            args.run_roms = False # never run ROMS if we skipped making the dot_in
        
        if args.run_roms:
            tt0 = time()
            # Run ROMS using the batch script.
            cmd_list = ['sbatch', '-p', args.cpu_choice, '-A', args.group_choice,
                str(roms_out_dir / 'klone_batch.sh')]
            proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # now we need code to wait until the run has completed
            
            # these are for checking on the run using squeue
            cmd_list = ['squeue', '-p', args.cpu_choice, '-A', args.group_choice]

            # first figure out if it has started
            for rrr in range(10):
                if rrr == 9:
                    print('Took too long for job to start: quitting')
                    sys.exit()
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                if jobname not in stdout.decode():
                    if args.verbose:
                        print('still waiting for run to start ' + str(rrr))
                        sys.stdout.flush()
                elif jobname in stdout.decode():
                    if args.verbose:
                        print('run started ' + str(rrr))
                        sys.stdout.flush()
                    break
                sleep(60)

            # and then figure out if it has finished (keeps looking for two hours, not always enough?)
            for rrr in range(60):
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                messages(stdout, stderr, 'Finished?', args.verbose)
                if jobname in stdout.decode():
                    if args.verbose:
                        print('still waiting ' + str(rrr))
                        sys.stdout.flush()
                    else:
                        pass
                elif (jobname not in stdout.decode()) and (len(stderr) == 0):
                    print(' - time to run ROMS = %d sec' % (time()-tt0))
                    sys.stdout.flush()
                    break
                else:
                    pass
                sleep(120)
    
            # The code here is a lot of claptrap to make sure there is a log file
            
            # A bit of checking to make sure that the log file exists...
            lcount = 0
            while not log_file.is_file():
                sleep(10)
                if args.verbose:
                    print(' - lcount = %d' % (lcount))
                    sys.stdout.flush()
                lcount += 1
                
                # trap for possible sbatch errors
                if lcount >= 10:
                    print(' - too long to write log.txt, assume there was some sbatch error')
                    sys.exit()
                    
            # ...and that it is done being written to.
            llcount = 0
            log_done = False
            while log_done == False:
                sleep(3)
                cmd_list = ['lsof', '-u', local_user,'|','grep',str(log_file)]
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                if args.verbose:
                    print(' - llcount = %d' % (llcount))
                    sys.stdout.flush()
                llcount += 1
                if str(log_file) not in stdout.decode():
                    if args.verbose:
                        print(' - log done and closed')
                        sys.stdout.flush()
                    log_done = True
                    
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
            for rrr in range(10):
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                if len(stderr) == 0: # it worked
                    break
                else:
                    sleep(20) # try again
            messages(stdout, stderr, 'Make output directory on ' + remote_machine, args.verbose)
            # (ii) move the contents of roms_out_dir
            cmd_list = ['scp','-r',str(roms_out_dir),
                remote_user + '@' + remote_machine + ':' + remote_dir0 + '/LO_roms/' + Ldir['gtagex']]
            for rrr in range(10):
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                if len(stderr) == 0: # it worked
                    break
                else:
                    sleep(20) # try again
            messages(stdout, stderr, 'Copy ROMS output to ' + remote_machine, args.verbose)
            # (iii) delete roms_out_dir and forcing files from several days in the past
            dt_prev = dt - timedelta(days=4)
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
    



