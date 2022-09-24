"""
This is the beginnings of a driver to run a bunch of tracker runs simultaneously.

It is meant to be a replacement for running a bunch of releases one after the other.

One issue I need to think through is the naming of things.

"""

import argparse
from lo_tools import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi

# command line arguments, can be input in any order
parser = argparse.ArgumentParser()

# Set the experiment name
# (details set in experiments.py, or, if it exists, user_experiments.py)

parser.add_argument('-gtx', '--gtagex', default='cas6_v0_live', type=str)
parser.add_argument('-ro', '--roms_out_num', default=0, type=int)
# 1 = Ldir['roms_out1'], etc.

# this is the first starting day
parser.add_argument('-d', '--date_string', default='2019.07.04', type=str)

parser.add_argument('-exp', '--exp_name', default='jdf0', type=str)
parser.add_argument('-clb', '--clobber', default=False, type=zfun.boolean_string)
# overwrite existing output folder if clobber == True
parser.add_argument('-sub_tag', default='', type=str)
# append an optional tag to the end of the output folder name

# These are False unless the flags are used with the argument True
# so if you do NOT use these flags the run will be:
# - trapped to the surface
# - no vertical turbulent diffusion
parser.add_argument('-3d', default=False, type=zfun.boolean_string) # do 3d tracking
parser.add_argument('-laminar', default=False, type=zfun.boolean_string) # no turbulence
parser.add_argument('-no_advection', default=False, type=zfun.boolean_string) # no advection
parser.add_argument('-sink', default=0, type=float) # particle sinking speed (m per day, e.g. 40)
parser.add_argument('-stay', default=0, type=float) # depth to try to stay at (m, e.g. 80)

# windage = a small number: 0 <= windage << 1 (e.g. 0.03)
# fraction of windspeed added to advection, only for 3d=False
parser.add_argument('-wnd', '--windage', default=0, type=float)

# You can make multiple releases using:
# number_of_start_days > 1 & days_between_starts, and which hour (UTC) to start on
parser.add_argument('-nsd', '--number_of_start_days', default=1, type=int)
parser.add_argument('-dbs', '--days_between_starts', default=1, type=int)
parser.add_argument('-dtt', '--days_to_track', default=1, type=int)
parser.add_argument('-sh', '--start_hour', default=0, type=int)

# number of divisions to make between saves for the integration
# e.g. if ndiv = 12 and we have hourly saves, we use a 300 sec step
# for the integration. 300 s seems like a good default value,
# based on Banas et al. (2009, CSR RISE paper).
parser.add_argument('-ndiv', default=12, type=int)
parser.add_argument('-sph', default=1, type=int)
# sph = saves per hour, a new argument to allow more frequent writing of output.

# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)

args = parser.parse_args()
argsd = args.__dict__

gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

# ==============================================================

# now run the releases in parallel

# make the list of start days (datetimes) for separate releases
idt_list = []
dt = datetime.strptime(TR['date_string'], '%Y.%m.%d')
for nic in range(TR['number_of_start_days']):
    idt_list.append(dt)
    dt = dt + timedelta(TR['days_between_starts'])


# do the initial extractions
N = args.number_of_start_days
proc_list = []
tt0 = time()
print('Working on ' + args.exp_name + ' (' + str(N) + ' start days)')
sys.stdout.flush()

for ii in range(N):
    cmd_list1 = ['python','tracker.py', 'and all the rest']
    # note that each job will hav 1 start day, and we need to increment these
    # in this driver
    proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
    proc_list.append(proc)

    # screen output about progress
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == N-1):
        print(str(ii))
        sys.stdout.flush()

    # Nproc controls how many subprocesses we allow to stack up
    # before we require them all to finish.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            proc.communicate()
        # make sure everyone is finished before continuing
        proc_list = []
    ii += 1
print(' Time to run all start days = %0.2f sec' % (time()- tt0))
sys.stdout.flush()
