"""
Shared helper functions for the forcing code, especially for argument passing.

"""
import argparse
import sys
from pathlib import Path
import Lfun, zfun # path to alpha set by calling function

def intro():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument('-g', '--gridname', type=str)   # e.g. cas6
    parser.add_argument('-t', '--tag', type=str)        # e.g. v3
    parser.add_argument('-f', '--frc', type=str)        # e.g. tide
    parser.add_argument('-r', '--run_type', type=str)   # backfill or forecast
    parser.add_argument('-s', '--start_type', type=str) # new or continuation
    parser.add_argument('-d', '--date_string', type=str) # e.g. 2019.07.04
    # optional arguments
    parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
    
    # get the args
    args = parser.parse_args()
    
    # test that required arguments were provided
    argsd = args.__dict__
    for a in ['gridname', 'tag', 'frc', 'run_type', 'start_type', 'date_string']:
        if argsd[a] == None:
            print('*** Missing required argument to forcing_argfun.intro(): ' + a)
            sys.exit()
        
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag)
    # add more entries to Ldir for use by make_forcing_main.py
    for a in ['frc', 'run_type', 'start_type', 'date_string', 'testing']:
        Ldir[a] = argsd[a]
        
    # create the expected output directories if needed
    # (a convenience when running make_forcing_main.py on its own while testing)
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']
    Lfun.make_dir(out_dir)
    Lfun.make_dir(out_dir / 'Info')
    Lfun.make_dir(out_dir / 'Data')

    return Ldir.copy()
    
def finale(Ldir, result_dict):
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']
    time_format = '%Y.%m.%d %H:%M:%S'
    s1 = '%s, %s\n' % ('start', result_dict['start_dt'].strftime(time_format))
    s2 = '%s, %s\n' % ('end', result_dict['end_dt'].strftime(time_format))
    s3 = '%s, %s\n' % ('result', result_dict['result'])
    if out_dir.is_dir():
        with open(out_dir / 'Info' / 'results.txt', 'w') as ffout:
            ffout.write(s1 + s2 + s3)
    else:
        print(' results.txt '.center(60,'+'))
        print(s1 + s2 + s3)
    

