"""
Shared helper functions for the dot_in code, especially for argument passing.

The new argument "tag_alt" is intended to facilitate different versions of a
ROMS run, where you change one piece of the forcing.  The dot_in code will
look to [gridname]_[tag] for forcing, but all other aspects will refer to
[gridname]_[tag_alt]_[ex_name], including the naming of the folder in LO/dot_in,
and the path to ROMS history files.

"""
import argparse
import sys
from pathlib import Path
import Lfun # path to alpha set by calling function

def intro():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument('-g', '--gridname', type=str)   # e.g. cas6
    parser.add_argument('-t', '--tag', type=str)        # e.g. v3
    parser.add_argument('-x', '--ex_name', type=str)    # e.g. lo8
    parser.add_argument('-r', '--run_type', type=str)   # backfill or forecast
    parser.add_argument('-s', '--start_type', type=str) # new or continuation
    parser.add_argument('-d', '--date_string', type=str) # e.g. 2019.07.04
    parser.add_argument('-bu', '--blow_ups', type=int) # e.g. 0
    parser.add_argument('-np', '--np_num', type=int) # e.g. 196, number of cores
    # optional arguments
    parser.add_argument('-ta', '--tag_alt', default='', type=str) # e.g. v3a
    
    # get the args
    args = parser.parse_args()
    
    # test that required arguments were provided
    argsd = args.__dict__
    for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'date_string',
                'blow_ups', 'np_num']:
        if argsd[a] == None:
            print('*** Missing required argument to forcing_argfun.intro(): ' + a)
            sys.exit()
        
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)
    
    # set tag_alt to tag if it is not provided
    if len(args.tag_alt) == 0:
        argsd['tag_alt'] = argsd['tag']
    Ldir['gtagex_alt'] = Ldir['gridname'] + '_' + argsd['tag_alt'] + '_' + Ldir['ex_name']
    
    # add more entries to Ldir for use by make_dot_in.py
    for a in ['run_type', 'start_type', 'date_string', 'blow_ups', 'np_num', 'tag_alt']:
        Ldir[a] = argsd[a]
        
    return Ldir.copy()
    
    

