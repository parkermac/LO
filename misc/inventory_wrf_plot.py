"""
Plot the inventory of WRF files.
"""

from datetime import datetime, timedelta
from lo_tools import Lfun
import pandas as pd
Ldir = Lfun.Lstart()
from pathlib import Path
from time import time

testing = False
if 'mac' in Ldir['lo_env']:
    wrf_dir = Ldir['data'] / 'wrf'
    testing = True
elif 'apogee' in Ldir['lo_env']:
    wrf_dir = Path('/pgdat1/parker/LO_data/wrf')
elif 'perigee' in Ldir['lo_env']:
    wrf_dir = Path('/data1/parker/LO_data/wrf')

if testing:
    year0 = 2019
    year1 = 2019
    dti = pd.date_range(start=datetime(year0,7,1), end=datetime(year1,7,7))
else:
    year0 = 2012
    year1 = datetime.now().year
    dti = pd.date_range(start=datetime(year0,1,1), end=datetime(year1,12,31))

tt0 = time()
df = pd.DataFrame(index=dti,columns=['d2','d3','d4'])
for dt in dti:
    fdir = wrf_dir / (dt.strftime('%Y%m%d')+'00')
    if fdir.is_dir():
        for dd in ['d2','d3','d4']:
            ddl = list(fdir.glob('*'+dd+'*'))
            ddl.sort()
            df.loc[dt,dd] = len(ddl)
    else:
        pass
print('Total time = %0.1f sec' % (time()-tt0))

out_dir = Ldir['LOo'] / 'misc'
Lfun.make_dir(out_dir)
out_fn = out_dir / 'inventory_wrf.p'
print('DataFrame saved to ' + str(out_fn))
df.to_pickle(out_fn)