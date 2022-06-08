"""
Shared helper functions for the dot_in code, especially for argument passing.

"""
import argparse
import sys
from pathlib import Path

pth = Path(__file__).absolute().parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

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
    parser.add_argument('-short_roms', default=False, type=Lfun.boolean_string)
    
    # get the args
    args = parser.parse_args()
    
    # test that required arguments were provided
    argsd = args.__dict__
    for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'date_string',
                'blow_ups', 'np_num']:
        if argsd[a] == None:
            print('*** Missing required argument to dot_in_argfun.intro(): ' + a)
            sys.exit()
        
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)
        
    # add more entries to Ldir for use by make_dot_in.py
    for a in ['run_type', 'start_type', 'date_string', 'blow_ups', 'np_num', 'short_roms']:
        Ldir[a] = argsd[a]
        
    return Ldir.copy()
    
    

