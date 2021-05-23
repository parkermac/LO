"""
This code moves/renames/copies files from old ROMS runs that were a large pile
of numbered history files for each save, and saves them in LO "day_folder" format.
"""

from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

import netCDF4 as nc
from datetime import datetime, timedelta

Ldir = Lfun.Lstart()

# case of the old ainlet files that Dave Sutherland made
in_dir = Ldir['parent'] / 'LiveOcean_roms' / 'output' / 'ainlet'

out_dir = Ldir['parent'] / 'LiveOcean_roms' / 'output' / 'ainlet_v0_n0'
Lfun.make_dir(out_dir)

fn_list = [ff for ff in in_dir.glob('ocean_his*nc')]
fn_list.sort()

dt_dict = dict()
for fn in fn_list:
    ds = nc.Dataset(fn)
    ot = ds['ocean_time']
    ots = ot[:].data[0]
    # assume units are like "seconds since 2006-01-01 00:00:00"
    otu = ot.units.split(' ')
    dt0 = datetime.strptime(otu[-2] + ' ' + otu[-1], '%Y-%m-%d %H:%M:%S')
    dt = dt0 + timedelta(seconds=ots)
    dt_dict[fn] = dt
    
# determine time between saves
delt = dt_dict[fn_list[1]] - dt_dict[fn_list[0]]
delth = delt.seconds/3600
fpd = 24/delth # files per day (not counting repeats)

fname_list = []
out_dict = dict()
for fn in fn_list:
    dt = dt_dict[fn]
    fname = 'f' + dt.strftime(Lfun.ds_fmt)
    if fname not in fname_list:
        fname_list.append(fname)
    hnum = (fpd/24)*(dt.hour + dt.minute/60) + 1
    hname = 'ocean_his_' + ('0000' + str(int(hnum)))[-4:] + '.nc'
    out_dict[fn] = out_dir / fname / hname
    
for fn in fn_list:
    print('%s => %s/%s' % (fn.name, out_dict[fn].parent.name, out_dict[fn].name))
    
    
# also need to copy 0001 to 0049 in previous day
