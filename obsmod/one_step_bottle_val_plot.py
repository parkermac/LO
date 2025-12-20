"""
Driver code to do a cast extraction, combine the results with bottle data, and plot the results to a png.

Testing on mac:
run one_step_bottle_val_plot -gtx cas7_t0_x4b -lt hourly -year0 2017 -test True
"""

import argparse
from time import time
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import sys

from lo_tools import Lfun

parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t1_x11ab
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
parser.add_argument('-lt', '--list_type', type=str, default='average') # hourly, average
parser.add_argument('-sources', type=str, default='all') # e.g. all, or other user-defined list
parser.add_argument('-otype', type=str, default='bottle') # observation type, e.g. ctd, bottle, etc.
parser.add_argument('-year0', type=int) # e.g. 2019
parser.add_argument('-year1', type=int, default=0) # will set to year0 if not specified
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex', 'roms_out_num', 'list_type', 'year0']:
    if argsd[a] == None:
        print('*** Missing required argument to extract_argfun.intro(): ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
# set year range
if Ldir['year1'] == 0:
    Ldir['year1'] = Ldir['year0'] + 1
year_list = list(range(Ldir['year0'],Ldir['year1']))

# Get the list of obs sources to use
if Ldir['sources'] == 'all':
    source_list = ['dfo1', 'ecology_nc', 'nceiSalish', 'nceiCoastal',
        'LineP', 'nceiPNW', 'WOD', 'kc', 'kc_pointJefferson']

# Loop over years:
for year in year_list:
    print('\n' + str(year) + '\n')

    # - Do the cast extractions for each source
    for source in source_list:
        cmd_list = ['python', str(Ldir['LO'] / 'extract' / 'cast' / 'extract_casts_fast.py'),
        '-gtx', Ldir['gtagex'],
        '-ro', str(Ldir['roms_out_num']),
        '-lt', Ldir['list_type'],
        '-source', source, 
        '-otype', Ldir['otype'],
        '-year', str(year)]

        if Ldir['testing']:
            print(cmd_list)
        else:
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            print(' ' + source)
            if len(stdout) > 0:
                a = stdout.decode()
                aa = a.split('\n')
                print('  %d casts processed' % (len(aa) - 1))
            else:
                print('  no casts found')
            if len(stderr) > 0:
                print('\n' + ' stderr '.center(60,'-'))
                print(stderr.decode())

# - Combine the cast extractions with obs values into a single DataFrame

# - Plot the results and save as a png.
