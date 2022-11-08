"""
Module of functions for LO version of LiveOcean.

This is purposely kept to a minimum of imports so that it will run with
whatever python3 exists on the large clusters we use for ROMS, e.g. mox and klone.

Re-coded 2022.04.19 to first look for LO_user/get_lo_info.py, and if not found
then look for LO/get_user_info.py.  The goal is to clean out LO_user of all of the code
that the LO developer (MacCready) would edit.  Then we add "hooks" to look for
user versions in LO_user at strategic places, as is done below.
"""
import os, sys, shutil
from pathlib import Path 
from datetime import datetime, timedelta
import importlib.util

# get initial version of Ldir when this module is loaded
upth = Path(__file__).absolute().parent.parent.parent.parent / 'LO_user'
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO'

# get the job_definitions module, looking first in LO_user
if (upth / 'get_lo_info.py').is_file():
    if str(upth) not in sys.path:
        sys.path.append(str(upth))
    try:
        import get_lo_info as glo
    except Exception as e:
        print(e)
elif (pth / 'get_lo_info.py').is_file():
    if str(pth) not in sys.path:
        sys.path.append(str(pth))
    try:
        import get_lo_info as glo
    except Exception as e:
        print(e)
else:
    print('Error from Lfun: missing LO/get_lo_info.py and LO_user/get_lo_info.py')
    sys.exit()

# initialize Ldir for this module
Ldir = glo.Ldir0.copy()

# this it the one place where the model time reference is set
modtime0 = datetime(1970,1,1,0,0)

# correct string for time units in ROMS forcing files
# see notes in Evernote, Run Log 9, 2020.10.06
roms_time_units = 'seconds since 1970-01-01 00:00:00'

# format used for naming day folders
ds_fmt = '%Y.%m.%d'

def Lstart(gridname='BLANK', tag='BLANK', ex_name='BLANK'):
    """
    This adds more run-specific entries to Ldir.
    """
    # put top level information from input into a dict
    Ldir['gridname'] = gridname
    Ldir['tag'] = tag
    Ldir['ex_name'] = ex_name
    # and add a few more things
    Ldir['gtag'] = gridname + '_' + tag
    Ldir['gtagex'] = gridname + '_' + tag + '_' + ex_name
    Ldir['grid'] = Ldir['data'] / 'grids' / gridname
    Ldir['forecast_days'] = 3
    Ldir['ds_fmt'] = ds_fmt
    Ldir['roms_time_units'] = roms_time_units
    Ldir['modtime0'] = modtime0
    return Ldir.copy()
    # the use of copy() means different calls to Lstart (e.g. when importing
    # plotting_functions) to not overwrite each other

def make_dir(pth, clean=False):
    """
    >>> WARNING: Be careful! This can delete whole directory trees. <<<
    
    Make a directory from the path "pth" which can:
    - be a string or a pathlib.Path object
    - be a relative path
    - have a trailing / or not
    
    Use clean=True to clobber the existing directory (the last one in pth).
    
    This function will create all required intermediate directories in pth.
    """
    if clean == True:
        shutil.rmtree(str(pth), ignore_errors=True)
    Path(pth).mkdir(parents=True, exist_ok=True)

def datetime_to_modtime(dt):
    """
    This is where we define how time will be treated
    in all the model forcing files.

    INPUT: dt is a single datetime value
    OUTPUT: dt as seconds since modtime0 (float)
    """
    t = (dt - modtime0).total_seconds()
    return t

def modtime_to_datetime(t):
    """
    INPUT: seconds since modtime0 (single number)
    OUTPUT: datetime version
    """
    dt = modtime0 + timedelta(seconds=t)
    return dt

def modtime_to_mdate_vec(mt_vec):
    """ 
    INPUT: numpy vector of seconds since modtime0
    - mt stands for model time
    OUTPUT: a vector of mdates
    """ 
    import matplotlib.dates as mdates
    #first make a list of datetimes
    dt_list = []
    for mt in mt_vec:
        dt_list.append(modtime0 + timedelta(seconds=mt))
    md_vec = mdates.date2num(dt_list)
    return md_vec
    
def boolean_string(s):
    # used by argparse (also present in zfun, redundant but useful)
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# Functions used by postprocessing code like pan_plot or the various extractors

def date_list_utility(dt0, dt1, daystep=1):
    """
    INPUT: start and end datetimes
    OUTPUT: list of LiveOcean formatted dates
    """
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime(ds_fmt))
        dt = dt + timedelta(days=daystep)
    return date_list

def fn_list_utility(dt0, dt1, Ldir, hourmax=24):
    """
    INPUT: start and end datetimes
    OUTPUT: list of all history files expected to span the dates
    - list items are Path objects
    """
    dir0 = Ldir['roms_out'] / Ldir['gtagex']
    fn_list = []
    date_list = date_list_utility(dt0, dt1)
    # new scheme 2022.10.09 to work with perfect restart
    dt00 = (dt0 - timedelta(days=1))
    fn_list.append(dir0 / ('f'+dt00.strftime(ds_fmt)) / 'ocean_his_0025.nc')
    for dl in date_list:
        f_string = 'f' + dl
        hourmin = 1
        for nhis in range(hourmin+1, hourmax+2):
            nhiss = ('0000' + str(nhis))[-4:]
            fn = dir0 / f_string / ('ocean_his_' + nhiss + '.nc')
            fn_list.append(fn)
    return fn_list
    
def get_fn_list(list_type, Ldir, ds0, ds1, his_num=2):
    """
    INPUT:
    A function for getting lists of history files.
    List items are Path objects
    """
    dt0 = datetime.strptime(ds0, ds_fmt)
    dt1 = datetime.strptime(ds1, ds_fmt)
    dir0 = Ldir['roms_out'] / Ldir['gtagex']
    if list_type == 'snapshot':
        # a single file name in a list
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = [dir0 / ('f' + ds0) / ('ocean_his_' + his_string + '.nc')]
    elif list_type == 'hourly':
        # list of hourly files over a date range
        fn_list = fn_list_utility(dt0,dt1,Ldir)
    elif list_type == 'daily':
        # list of history file 21 (Noon PST) over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'weekly':
        # like "daily" but at 7-day intervals
        fn_list = []
        date_list = date_list_utility(dt0, dt1, daystep=7)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'allhours':
        # a list of all the history files in a directory
        # (this is the only list_type that actually finds files)
        in_dir = dir0 / ('f' + ds0)
        fn_list = [ff for ff in in_dir.glob('ocean_his*nc')]
        fn_list.sort()

    return fn_list
    
def choose_item(in_dir, tag='', exclude_tag='',
    itext='** Choose item from list **', last=False):
    """
    INPUT: in_dir = Path object of a directory
    OUTPUT: just the name you chose, a string, not the full path.
    
    You can set strings to search for (tag), strings to exclude (exclude_tag),
    and the prompt text (itext).
    
    Use last=True to have the "return" choice be the last one.
    """
    print('\n%s\n' % (itext))
    ilist_raw = [item.name for item in in_dir.glob('*')]
    ilist_raw.sort()
    if len(tag) == 0:
        ilist = [item for item in ilist_raw if item[0] != '.']
    else:
        ilist = [item for item in ilist_raw if tag in item]
        
    if len(exclude_tag) == 0:
        pass
    else:
        ilist = [item for item in ilist if exclude_tag not in item]
    
    Nitem = len(ilist)
    idict = dict(zip(range(Nitem), ilist))
    for ii in range(Nitem):
        print(str(ii) + ': ' + ilist[ii])
    if last == False:
        my_choice = input('-- Input number -- (return=0) ')
        if len(my_choice)==0:
            my_choice = 0
    elif last == True:
        my_choice = input('-- Input number -- (return=last) ')
        if len(my_choice)==0:
            my_choice = Nitem-1
        
    my_item = idict[int(my_choice)]
    return my_item
    
def dict_to_csv(in_dict, out_fn):
    """
    Writes a dict to a csv file.
    
    out_fn should be a Path object.
    """
    out_fn.unlink(missing_ok=True)
    with open(out_fn, 'w') as f:
        for k in in_dict.keys():
            f.write(k + ',' + str(in_dict[k]) + '\n')
            
def csv_to_dict(in_fn):
    """
    Reads a csv file into a dict and returns it.
    
    We should add some error checking to make sure the input is as expected.
    """
    out_dict = dict()
    with open(in_fn, 'r') as f:
        for line in f:
            k,v = line.split(',')
            out_dict[k] = str(v).replace('\n','')
    return out_dict
    
def messages(stdout, stderr, mtitle, test_flag=True):
    """
    A utility function for displaying subprocess info.
    """
    if test_flag:
        print((' ' + mtitle + ' ').center(60,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print(' stderr '.center(60,'-'))
            print(stderr.decode())
        sys.stdout.flush()
        
def module_from_file(module_name, file_path):
    """
    This is used for the hook ot LO_user.  It allows you to import a module from a
    specific path, even if a module of the same name exists in the current directory.
    """
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

if __name__ == '__main__':
    # TESTING: run Lfun will execute these
    
    if False:
        print(' TESTING Lstart() '.center(60,'-'))
        Ldir = Lstart(gridname='cas6', tag='v3', ex_name='lo8b')
        print(' Ldir seen by make_forcing_main '.center(60,'+'))
        for k in Ldir.keys():
            print('%20s : %s' % (k, Ldir[k]))
            
    if False:
        print(' TESTING datetime_to_modtime() '.center(60,'-'))
        dt = datetime(2019,7,4)
        print(dt)
        t = datetime_to_modtime(dt)
        print(t)
        print(' TESTING modtime_to_datetime() '.center(60,'-'))
        dt_new = modtime_to_datetime(t)
        print(dt_new)
    
    if False:
        print(' TESTING copy_to_azure() '.center(60,'-'))
        input_filename = Ldir['data'] / 'accounts' / 'az_testfile.txt'
        output_filename = input_filename.name
        container_name = 'pm-share'
        az_dict = copy_to_azure(input_filename, output_filename, container_name, Ldir)
        if az_dict['result'] =='success':
            print('USE THIS URL TO ACCESS THE FILE')
            print(az_dict['az_url'])
        elif az_dict['result'] =='fail':
            print('EXCEPTION')
            print(az_dict['exception'])
            
    if True:
        print(' TESTING get_fn_list() '.center(60,'-'))
        Ldir = Lstart(gridname='cas6', tag='v0', ex_name='live')
        # Ldir['roms_out'] = Ldir['roms_out']
        list_type = 'allhours'
        ds0 = '2019.07.04'
        ds1 = '2019.07.05'
        for list_type in ['daily','snapshot', 'allhours', 'hourly']:
            print(list_type.center(60,'.'))
            fn_list = get_fn_list(list_type, Ldir, ds0, ds1, his_num=7)
            for fn in fn_list:
                print(fn)
                
    if False:
        print(' TESTING choose_item() '.center(60,'-'))
        Ldir = Lstart(gridname='cas6', tag='v3', ex_name='lo8b')
        in_dir = Ldir['roms_out1'] / Ldir['gtagex'] / 'f2019.07.04'
        my_item = choose_item(in_dir, tag='.nc', exclude_tag='')
        print(my_item)
    
    



