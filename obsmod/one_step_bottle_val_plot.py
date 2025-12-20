"""
Driver code to do a cast extraction, combine the results with bottle data, and plot the results to a png.
"""

import argparse
from lo_tools import Lfun

parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t1_x11ab
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
parser.add_argument('-lt', '--list_type', type=str, default='average') # hourly, average
parser.add_argument('-source_list', type=str, default='all') # e.g. all, or other user-defined list
parser.add_argument('-otype', type=str, default='bottle') # observation type, e.g. ctd, bottle, etc.
parser.add_argument('-year0', type=int) # e.g. 2019
parser.add_argument('-year1', type=int, default=0) # will set to year0 if not specified
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex']:
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

# Get the list of obs sources to use

# Loop over years:

# - Do the cast extractions for each source

# - Combine the cast extractions with obs values into a single DataFrame

# - Plot the results and save as a png.
