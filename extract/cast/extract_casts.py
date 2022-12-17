"""
This is code for doing cast extractions.

Refactored 2022_07 to conform to the new cast data format.

Test on mac in ipython:
run extract_casts -gtx cas6_v0_live -source dfo -otype ctd -year 2019 -test True

"""

import sys
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime

from lo_tools import Lfun, zfun, zrfun
from lo_tools import extract_argfun as exfun
import cast_functions as cfun

Ldir = exfun.intro() # this handles the argument passing

year_str = str(Ldir['year'])

out_dir = (Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'cast' /
    (Ldir['source'] + '_' + Ldir['otype'] + '_' + year_str))
Lfun.make_dir(out_dir, clean=True)

info_fn = Ldir['LOo'] / 'obs' / Ldir['source'] / Ldir['otype'] / ('info_' + year_str + '.p')
if info_fn.is_file():
    ii = 0
    info_df = pd.read_pickle(info_fn)
    for cid in info_df.index:
        lon = info_df.loc[cid,'lon']
        lat = info_df.loc[cid,'lat']
        dt = info_df.loc[cid,'time']
        out_fn = out_dir / (str(int(cid)) + '.nc')
        fn = cfun.get_his_fn_from_dt(Ldir, dt)
                
        if fn.is_file(): # useful for testing
            
            # check on which bio variables to get
            if ii == 0:
                ds = xr.open_dataset(fn)
                if 'NH4' in ds.data_vars:
                    npzd = 'new'
                elif 'NO3' in ds.data_vars:
                    npzd = 'old'
                else:
                    npzd = 'none'
                ds.close()
            
            print('Get ' + out_fn.name)
            sys.stdout.flush()
            cfun.get_cast(out_fn, fn, lon, lat, npzd)
            ii += 1
            if Ldir['testing'] and (ii > 3):
                break
        else:
            pass



