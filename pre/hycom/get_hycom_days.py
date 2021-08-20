"""
Extract and save extracted fields from a sequence of hycom past days.

"""

from lo_tools import Lfun, zfun
from lo_tools import hycom_functions as hfun

import pickle
import netCDF4 as nc
from datetime import datetime, timedelta

Ldir = Lfun.Lstart()

# optional command line input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-force', '--force_overwrite', default=False, type=zfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
args = parser.parse_args()
force_overwrite = args.force_overwrite
testing = args.testing

# initial experiment list
h_list = list(hfun.hy_dict.keys())
h_list.sort()

if testing:
    h_list = h_list[-2:]
    
if testing:
    var_list = 'surf_el'
else:
    var_list = 'surf_el,water_temp,salinity,water_u,water_v'

# specify output directory
out_dir0 = Ldir['data'] / 'hycom'
dt_dir = out_dir0 / 'dt_lists'

f = open(out_dir0 / 'log.txt', 'w+')

# loop over all days in all lists
for hy in h_list:
    print('\n** Working on ' + hy + ' **')
    f.write('\n\n** Working on ' + hy + ' **')
    
    out_dir = out_dir0 / hy
    Lfun.make_dir(out_dir)
    
    dt_list = pickle.load(open(dt_dir / (hy + '.p'), 'rb'))
    
    if testing:
        dt_list = dt_list[:2]
        
    for dt in dt_list:
        print(' - ' + datetime.strftime(dt, Lfun.ds_fmt))
        
        out_fn = out_dir / ('h' + datetime.strftime(dt, Lfun.ds_fmt) + '.nc')
        
        if force_overwrite:
            out_fn.unlink(missing_ok=True)
                
        if not out_fn.is_file():
            result = hfun.get_extraction(hy, dt, out_fn, var_list)
            f.write('\n ' + datetime.strftime(dt, Lfun.ds_fmt) + ' ' + result)
            
f.close()