"""
Shared helper functions for the forcing code, especially for argument passing.

"""
import argparse
import sys
from lo_tools import Lfun

def intro():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument('-g', '--gridname', type=str)   # e.g. cas6
    parser.add_argument('-f', '--frc', type=str)        # e.g. tide
    parser.add_argument('-r', '--run_type', type=str)   # backfill or forecast
    parser.add_argument('-s', '--start_type', type=str, default='continuation') # new, continuation, or perfect
    parser.add_argument('-d', '--date_string', type=str) # e.g. 2019.07.04
    # optional arguments
    parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
    
    # optional arguments used only for ocnN, to determine what to nest inside
    parser.add_argument('-gtx', '--gtagex', default='cas6_traps2_x2b', type=str) # e.g. cas6_traps2_x2b
    parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
    parser.add_argument('-do_bio', default=False, type=Lfun.boolean_string) # True to add bio vars to ocnN forcing
    
    # get the args
    args = parser.parse_args()
    
    # test that required arguments were provided
    argsd = args.__dict__
    for a in ['gridname', 'frc', 'run_type', 'start_type', 'date_string']:
        if argsd[a] == None:
            print('*** Missing required argument to forcing_argfun.intro(): ' + a)
            sys.exit()
        
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=args.gridname)
    # add more entries to Ldir for use by make_forcing_main.py
    for a in ['frc', 'run_type', 'start_type', 'date_string', 'testing','gtagex','roms_out_num','do_bio']:
        Ldir[a] = argsd[a]
    # set where to look for model output
    if Ldir['roms_out_num'] == 0:
        pass
    elif Ldir['roms_out_num'] > 0:
        Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
        
    # create the expected output directories if needed
    # (a convenience when running make_forcing_main.py on its own while testing)
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']
    Lfun.make_dir(out_dir)
    Lfun.make_dir(out_dir / 'Info')
    Lfun.make_dir(out_dir / 'Data')

    return Ldir.copy()
    
def finale(Ldir, result_dict):
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']
    time_format = '%Y.%m.%d %H:%M:%S'
    total_sec = (result_dict['end_dt']-result_dict['start_dt']).total_seconds()
    if 'note' in result_dict.keys():
        pass
    else:
        result_dict['note'] = 'NONE'
        
    s1 = ('* frc=%s, day=%s, result=%s, note=%s\n' %
        (Ldir['frc'], Ldir['date_string'], result_dict['result'], result_dict['note']))
    
    s2 = ('  start=%s (took %d sec)\n' %
        (result_dict['start_dt'].strftime(time_format), int(total_sec)))
    
    s3 = ('  %s\n' % (str(out_dir)))
    
    with open(out_dir / 'Info' / 'results.txt', 'w') as ffout:
        ffout.write(s1 + s2 + s3)
    

