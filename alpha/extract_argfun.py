"""
Shared helper functions for the extraction code, especially for argument passing.

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
    parser.add_argument('-x', '--ex_name', type=str)    # e.g. lo8b
    # optional arguments
    parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = roms1, etc.
    parser.add_argument('-0', '--ds0', type=str, default='')        # e.g. 2019.07.04
    parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
    parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
    # optional extra arguments specific to different types of extractions
    parser.add_argument('-a1', type=str, default='')
    parser.add_argument('-a2', type=str, default='')
    parser.add_argument('-a3', type=str, default='')
    
    # get the args
    args = parser.parse_args()
    
    # test that required arguments were provided
    argsd = args.__dict__
    for a in ['gridname', 'tag', 'ex_name']:
        if argsd[a] == None:
            print('*** Missing required argument to forcing_argfun.intro(): ' + a)
            sys.exit()
        
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)
    # add more entries to Ldir for possible use by extractors
    for a in ['roms_out_num', 'ds0', 'ds1', 'testing', 'a1', 'a2', 'a3']:
        Ldir[a] = argsd[a]
        
    # set the end day string
    if len(Ldir['ds1'])==0:
        Ldir['ds1'] = Ldir['ds0']
        
    # set where to look for model output
    if args.roms_out_num == 0:
        pass
    elif args.roms_out_num > 0:
        Ldir['roms_out'] = Ldir['roms_out' + str(args.roms_out_num)]
        
    return Ldir.copy()
    

