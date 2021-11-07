"""
This is the main program for cleaning out unneeded history files from past forecasts. It
deletes history files 26-73 (or more if the forecast is longer) for the DAY BEFORE the one
it is passed as date_string.

Test on apogee:
python forecast_cleaner.py -gtx cas6_v0_u0kb -ro 0 -r backfill -0 [today's datestring] -test True
Testing just prints what it would do, but does not actually delete the files.

Test on mac:
run forecast_cleaner.py -gtx cas6_v3_lo8b -ro 2 -r backfill -0 2019.07.05 -test True

Run for real on apogee:
python forecast_cleaner.py -gtx cas6_v0_u0kb -ro 0 -r forecast > clean.log &

And of course it can be run for a single forecast day (e.g. in cron) or any past range (by hand).
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta
from lo_tools import Lfun
import argparse
from lo_tools import Lfun

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-r', '--run_type', type=str)   # forecast or backfill
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
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

# set time range to process
if args.run_type == 'forecast':
    ds0 = datetime.now().strftime(Lfun.ds_fmt)
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = dt0
    ds1 = dt1.strftime(Lfun.ds_fmt)
elif args.run_type == 'backfill':
    ds0 = args.ds0
    if len(args.ds1) == 0:
        ds1 = ds0
    else:
        ds1 = args.ds1
    dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
    dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
else:
    print('Error: Unknown run_type')
    sys.exit()
    
keep_name_list = []
for ii in range(1,26):
    iistr = ('0000' + str(ii))[-4:]
    keep_name_list.append('ocean_his_' + iistr + '.nc')

# loop over all days
print(' forecast_cleaner '.center(60,'='))
dt = dt0
while dt <= dt1:

    # the day before date_string
    dt_prev = dt - timedelta(days=1)
    ds_prev = dt_prev.strftime(Lfun.ds_fmt)

    f_string_prev = 'f' + ds_prev
    dir_prev = Ldir['roms_out'] / Ldir['gtagex'] / f_string_prev
    if dir_prev.is_dir():
        print('Deleting forecast files for ' + ds_prev)
        do_delete = True
    else:
        print('Directory not found for: ' + ds_prev)
        do_delete = False
    
    if do_delete:
        fn_list = [ff for ff in dir_prev.glob('ocean_his*nc')]
        fn_list.sort()
        jj = 0
        for fn in fn_list:
            if fn.name in keep_name_list:
                pass
            else:
                if Ldir['testing']:
                    print(' - would delete ' + str(fn))
                else:
                    fn.unlink(missing_ok=True)
                    jj += 1
        print(' - ' + str(jj) + ' files deleted')
        
    dt += timedelta(days=1)
    

