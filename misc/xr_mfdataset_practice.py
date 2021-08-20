"""
Testing the open_mfdataset function in xarray.

RESULT: This achieves my initial goal of opening a multi-file dataset,
and creating a weighted mean.  It works great for the 25 files in a day folder
but fails for the 73 files in a three-day sequence, even when I specify smaller chunks.
So at this point it is sort of useful, but not for proper tidal averaging.
"""

from lo_tools import Lfun, zrfun, zfun

import sys
import numpy as np
import pandas as pd
import xarray as xr
from time import time
import matplotlib.pyplot as plt

# open a ROMS history file
Ldir = Lfun.Lstart(gridname='cas6',tag='v3',ex_name='lo8b')
Ldir['roms_out'] = Ldir['roms_out2']
fn_list = Lfun.get_fn_list('hourly',Ldir,'2019.07.04','2019.07.04')

tt0 = time()
ds = xr.open_mfdataset(fn_list, concat_dim='ocean_time', parallel=True)#, chunks={'ocean_time':1, 's_rho':1, 's_w':1},
# chunks does not appear to speed things up in this case
print('Time to open mfdataset = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

NT = len(ds.ocean_time)
f = zfun.hanning_shape(NT)
# note that the filter does not have to sum to 1.  It appears that the xarray "weighted"
# method enforces that internally
tt0 = time()
ds['filt'] = (('ocean_time'), f)
print('Time to add filt = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
    
# make full-volume weighted averages
for vn in ['salt', 'temp', 'oxygen']:
    tt0 = time()
    m = ds[vn].weighted(ds.filt).mean().values
    print('%s: mean = %0.3f (%0.2f sec)'% (vn, m, (time()- tt0)))
    sys.stdout.flush()

# make the weighted average surface salinity
mf = ds.salt[:,-1,:,:].weighted(ds.filt).mean(dim='ocean_time')

plt.close('all')
mf.plot()
plt.show()
