"""
Extract fields at a number of sections which may be used later for TEF analysis
of transport and transport-weighted properties.

All input parameters specified at the command line, so this can be run in the background
because it can take a few hours.

Takes about 10-15 hours for 39 cas6 sections, per year.

To test on mac (default is to just get ai1 section):

run extract_sections -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -get_bio True -test True

To get all sections for the same time:

python extract_sections.py -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -sect_name all > log_all_test &

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun

Ldir = exfun.intro() # this handles the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# set list of variables to extract
if Ldir['get_bio']:
    vn_list = ['salt', 'temp', 'oxygen', 'NO3', 'TIC', 'alkalinity']
else:
    vn_list = ['salt']

import Lfun
import numpy as np
import netCDF4 as nc

import zrfun
import tef_fun

if Ldir['testing']:
    from importlib import reload
    reload(tef_fun)

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']
dt0 = datetime.strptime(ds0, Ldir['ds_fmt'])
dt1 = datetime.strptime(ds1, Ldir['ds_fmt'])
ndays = (dt1-dt0).days + 1

print('Working on:')
outname = 'extractions_' + ds0 + '_' + ds1
print(outname)

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / outname
Lfun.make_dir(out_dir, clean=True)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
# initialize a dictionary of info for each section
sect_info = dict()
# select which sections to extract
if Ldir['sect_name'] == 'all':
    # full list
    sect_list = [item for item in sect_df.index]
else: # single item
    if Ldir['sect_name'] in sect_df.index:
        sect_list = [Ldir['sect_name']]
    else:
        print('That section is not available')
        sys.exit()

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
NT = len(fn_list)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
NZ = S['N']

print('\nGetting section definitions and indices:')
for sect_name in sect_list:
    print(sect_name)
    # name output file
    out_fn = out_dir / (sect_name + '.nc')
    # get section lat, lon, and other info
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    # get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    NX = len(Mask) # this the length of the section, not specific to x or y
    # save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, NT, NX, NZ, out_fn)
    # initialize a netcdf file for this section
    tef_fun.start_netcdf(fn, out_fn, NT, NX, NZ, Lon, Lat, Ldir, vn_list)
    # note that this function deletes the existing out_fn

# extract and save time-dependent fields
count = 0
print('\nStarting extraction of fields:')
print(vn_list)
for fn in fn_list:
    if np.mod(count,24)==0:
        print('  working on %d of %d' % (count, NT))
        sys.stdout.flush()
    ds = nc.Dataset(fn)
    # loop over all sections
    for sect_name in sect_list:
        sinfo = sect_info[sect_name]
        # this is where we add the data from this history file
        # to all of the sections, each defined by sinfo
        tef_fun.add_fields(ds, count, vn_list, G, S, sinfo)
    ds.close()
    count += 1
    
# test for success 
if True: # placeholder for a test
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()



