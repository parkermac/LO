"""
Shared helper functions for the extraction code, especially for argument passing.

"""
import argparse
import sys
from lo_tools import Lfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from datetime import datetime

def intro():
    parser = argparse.ArgumentParser()
    # which run to use
    parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
    parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.
    # select time period and frequency
    parser.add_argument('-r', '--run_type', type=str, default='')   # backfill or forecast
    parser.add_argument('-d', '--date_string', type=str) # e.g. 2019.07.04
    parser.add_argument('-job', default='', type=str) # e.g. surface0
    parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
    parser.add_argument('-Nproc', type=int, default=10) # number of subprocesses to use
    # get the args and put into Ldir
    args = parser.parse_args()
    argsd = args.__dict__
    for a in ['gtagex', 'roms_out_num', 'date_string', 'job']:
        if argsd[a] == None:
            print('*** Missing required argument to post_argfun.intro(): ' + a)
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

    # create the expected output directories if needed
    # (a convenience when running post_main.py on its own while testing)
    out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
    Lfun.make_dir(out_dir)
    Lfun.make_dir(out_dir / 'Info')
    Lfun.make_dir(out_dir / 'Data')

    return Ldir.copy()
    
def finale(Ldir, result_dict):
    out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
    time_format = '%Y.%m.%d %H:%M:%S'
    total_sec = (result_dict['end_dt']-result_dict['start_dt']).total_seconds()
    
    s1 = ('%s %s %s\n' %
        (Ldir['job'], Ldir['date_string'], result_dict['result'].upper()))
    
    s2 = ('start=%s (took %d sec)\n' %
        (result_dict['start_dt'].strftime(time_format), int(total_sec)))
    
    s3 = ('%s\n' % (str(out_dir)))
    
    S = s1 + s2 + s3
    
    if 'note' in result_dict.keys():
        S = S + ('NOTE: %s\n' % (result_dict['note']))
    
    with open(out_dir / 'Info' / 'results.txt', 'w') as ffout:
        ffout.write(S)
        
def copy_to_server(Ldir, out_fn, subdir=''):
    """
    Copy the  extraction file to the server and write a little "done" file.
    
    The optional "subdir" argument allows writing to a subdirectory. The first use of this
    is for hourly files made by splitting up the harcourt extraction.
    """
    if ('apogee' in Ldir['lo_env']) and (Ldir['testing'] == False):
    
        share_user = 'parker@liveocean.apl.uw.edu'
        if len(subdir) == 0:
            share_dir = '/data/www/liveocean/output/' + 'f' + Ldir['date_string']
        else:
            share_dir = '/data/www/liveocean/output/' + 'f' + Ldir['date_string'] + '/' + subdir

        # (i) make the output directory
        cmd_list = ['ssh', share_user, 'mkdir -p ' + share_dir]
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        Lfun.messages(stdout, stderr, 'Make output directory on server')
    
        is_done = False
        try:
            # (ii) copy the extraction to there
            cmd_list = ['scp',str(out_fn), share_user + ':' + share_dir]
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            Lfun.messages(stdout, stderr, 'Copying extraction to ' + share_dir)
            is_done = True
        except Exception as e:
            # I wonder if checking stderr might be a better error trap?
            print('Problem moving file to server')
            print(e)
        
        if is_done and (len(subdir) == 0):
            # (iii) then write a little text file to alert users
            share_name = out_fn.name.replace('.nc','')
            out_dir = out_fn.parent
            done_fn = out_dir / (share_name + '_done.txt')
            done_fn.unlink(missing_ok=True)
            with open(done_fn, 'w') as ffout:
                ffout.write(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
            cmd_list = ['scp',str(done_fn), share_user + ':' + share_dir]
            proc = Po(cmd_list, stdout=Pi, stderr=Pi)
            stdout, stderr = proc.communicate()
            Lfun.messages(stdout, stderr, 'Copying done file to server')
            
    else:
        print('** Did not copy output to server! **')
        print('Warning: copying extractions to server for sharing only works for parker from apogee.')
        print('Also it will not copy when testing is True.')
        
