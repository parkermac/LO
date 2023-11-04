"""
Code to generate an inventory of WRF files.
"""

from datetime import datetime, timedelta
from lo_tools import Lfun
import pandas as pd
Ldir = Lfun.Lstart()
import matplotlib.pyplot as plt
from pathlib import Path

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

plt.close('all')
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(111)
df.plot(ax=ax)
fig.savefig('inventory_wrf.png')
