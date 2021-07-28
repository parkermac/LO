"""
Generic code to plot any mooring extraction - meant to function outside of the
LO system so it is easier for others to use.  Put the mooring file in the
same directory as this code.

"""

from pathlib import Path
from datetime import datetime, timedelta
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# this it the one place where the model time reference is set
modtime0 = datetime(1970,1,1,0,0)
def modtime_to_datetime(t):
    """
    INPUT: seconds since modtime0 (single number)
    OUTPUT: datetime version
    """
    dt = modtime0 + timedelta(seconds=t)
    return dt

# choose the file
in_dir = Path(__file__).absolute().parent
moor_fn = in_dir / 'WestPoint_2020.01.01_2020.12.31.nc'

# open the dataset
ds = nc.Dataset(moor_fn)

# get time axis
ot = ds['ocean_time'][:]
print('time step of mooring'.center(60,'-'))
print(np.diff(ot))
# and convert to a list of datetimes
tind = [modtime_to_datetime(tt) for tt in ot]
print('time limits'.center(60,'-'))
print('start ' + str(tind[0]))
print('end   ' + str(tind[-1]))

# print info about the variables
print('info'.center(60,'-'))
VN_list = []
for vn in ds.variables:
    try:
        print('%s (%s: %s) %s' % (vn, ds[vn].long_name, ds[vn].units, ds[vn].shape))
    except AttributeError:
        print('%s %s' % (vn, ds[vn].shape))
    VN_list.append(vn)
    
# populate list of variables to plot
vn_list = []
if 'salt' in VN_list:
    vn_list += ['salt', 'temp']
if 'NO3' in VN_list:
    vn_list += ['NO3', 'oxygen', 'phytoplankton']
if 'u' in VN_list:
    vn_list += ['u', 'v']

# put time series of variables at surface in a pandas DataFrame
df = pd.DataFrame(index=tind)
df['zeta'] = ds['zeta'][:,0,0]
# add surface variables
for vn in vn_list:
    df[vn] = ds[vn][:, -1, 0, 0]

# and plot it
plt.close('all')
df.plot(subplots=True, figsize=(16,10))
plt.show()

