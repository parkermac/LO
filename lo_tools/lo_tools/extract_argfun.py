"""
Shared helper functions for the extraction code, especially for argument passing.

"""
import argparse
import sys
from lo_tools import Lfun

def intro():
    parser = argparse.ArgumentParser()
    # which run to use
    parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v0_live
    parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
    # select time period and frequency
    parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
    parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
    parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
    parser.add_argument('-Nproc', type=int, default=10) # number of subprocesses to use
    # arguments used by extract/cast
    parser.add_argument('-source', type=str) # e.g. dfo
    parser.add_argument('-otype', type=str) # observation type, e.g. ctd, bottle, etc.
    parser.add_argument('-year', type=int) # e.g. 2019
    # arguments used by extract/tef and tef2
    parser.add_argument('-sect_name', type=str, default='ai1')
    parser.add_argument('-get_bio', type=Lfun.boolean_string, default=False)
    # arguments used by extract/tef2
    parser.add_argument('-ctag','--collection_tag', type=str)
    parser.add_argument('-riv', type=str) # e.g. riv00
    
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

    return Ldir.copy()
    

