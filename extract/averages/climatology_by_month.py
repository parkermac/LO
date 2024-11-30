"""
Code to create monthly climatology files by averaging the month_mean files
for a given month over a range of years.

This is hard-coded instead of being an LO command-line tool because:
(i) This is a task that will only be done occasionally.

This is quick to run: took about 15 minutes to do means of 10 years on apogee.
"""

import xarray as xr
from datetime import datetime, timedelta
from time import time
from pathlib import Path
import sys

from lo_tools import Lfun, zrfun

Ldir = Lfun.Lstart()

testing = False

gtagex = 'cas7_t0_x4b'
# For this run all the monthly means are in /dat1/parker/LO_roms/cas7_t0_x4b/averages
# with names like "month_mean_2014_01.nc".

years = range(2014,2024) # 2014-2023 is 10 years
averaging_factor = 1/len(years)

in_dir = Path('/dat1/parker/LO_roms/cas7_t0_x4b/averages')
out_dir = Path('/dat1/parker/LO_roms/cas7_t0_x4b/climatologies')
if testing == False:
    Lfun.make_dir(out_dir)

for month in range(1,13):
    tt0 = time()
    ii = 0
    for year in years:
        dt = datetime(year,month,1)
        this_ym = dt.strftime('%Y_%m')
        in_fn = in_dir / ('monthly_mean_' + this_ym + '.nc')
        if testing == False:
            ds = xr.open_dataset(in_fn, decode_times=False)
            if ii == 0:
                lp = (ds*averaging_factor).compute()
            else:
                lp = (lp + ds*averaging_factor).compute()
            ds.close()
            ii += 1
    # save the temp file to NetCDF
    fstr = ('000' + str(month))[-2:]
    out_fn = out_dir / ('monthly_clim_' + fstr + '.nc')
    if testing == False:
        # and save to NetCDF, with compression
        Enc_dict = {vn:zrfun.enc_dict for vn in lp.data_vars}
        lp.to_netcdf(out_fn, encoding=Enc_dict)
        # NOTE: for some reason this compression works fine in this case but
        # throws an error "OverflowError: Python int too large to convert to C long"
        # when I try to do it in extract_lowpass.py. No clue why.

    print('Month = ' + str(month))
    print(str(out_fn))
    print(' - Time to make monthly_clim = %0.1f minutes' % ((time()-tt0)/60))
    sys.stdout.flush()

