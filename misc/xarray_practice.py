"""
Code to practice and test the use of xarray in the LO system.
"""
    
from lo_tools import Lfun, zrfun

import numpy as np
import pandas as pd
import xarray as xr

# open a ROMS history file
Ldir = Lfun.Lstart()
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2019.07.04' / 'ocean_his_0001.nc')
G,S,T = zrfun.get_basic_info(fn)

ds = xr.open_dataset(fn)
"""
An xarray Dataset with properties like:
- dims (like xi_rho)
- coords (lke lon_rho)
- attrs
- data_vars (dict-like of all the variables, e.g. salt)

xr.open_dataset() = lazy loading
xr.load_dataset() = load everything

For a typical LiveOcean history file from cas6_v3_lo8b the original file is
2 GB, and loading it using load_dataset() appears to use about 700 MB of memory,
whereas open_dataset() only appears to use 1 MB!

Here are typical lines for making a Dataset from scratch and filling it:
ds = xr.Dataset(coords={'time': times,'seg': segs})
for vn in vn_list:
    v = some appropriate ndarray
    ds[vn] = (('time','seg'), v)
"""

s = ds.salt
"""
An xarray DataArray of shape (1, 30, 1302, 663)
with properties like:
- values (typically an ndarray)
- dims (like xi_rho)
- coords (lke lon_rho)
- attrs

NOTE: the coords are shared across all variables in the Dataset
"""

# NOTE: hereafter we refer to a DataArray as simply an array

# make a DataArray from scratch:
a = xr.DataArray(np.arange(12).reshape((3,4)), dims=['x','y'], coords={'x':[1,2,3], 'y':list('abcd')})
# There are many ways to make a DataArray.  In the version above we passed it data (required),
# coordinates as a dict, and dims.  Passing dims is required when specifying the coordinates
# as a dict.
# The version below is a different way of making the same array using a list of tuples for the
# coordinates
a_alt = xr.DataArray(np.arange(12).reshape((3,4)), [('x',[1,2,3]), ('y',list('abcd'))])
# One could also just pass a pandas Series for a 1D array.

# ways to work with it (these generally return another DataArray)
a = a * 2 # math
am = a.mean(dim='y')

# There are FOUR ways to index the array.
# - these two ways CAN change values in the array:
a[0, :] = 77 # positional and by integer label, like numpy
a.loc[2,'b'] = 88 # loc or "location": positional and coordinate label, like pandas
# - whereas these two ways can only return a new array, not change values in the array
aa = a.isel(x=0) # isel or "integer select":  by dimension name and integer label
bb = a.sel(y='b') # sel or "select": by dimension name and coordinate label, or a list of labels

# plotting
t = pd.date_range('1/1/2021','1/10/2021', freq='H') # xarray loves pandas time Index objects
# t = pd.date_range('2021-1-1','2021-1-10', freq='H') # equivalent
x = np.linspace(0,10,len(t))
A = xr.DataArray(np.sin(x), dims=('t'), coords={'t':t})
A.attrs['long_name'] = 'Things'
A.attrs['units'] = 'cubits'
# then A.plot() will make a nicely labeled line plot

# and here is where it gets really clever: broadcasting!
B = xr.DataArray([1,2,3], dims=('q'), coords={'q':list('LMN')})
# then these operations make 2D arrays from 1D:
C = A + B # the order of the dimensions appears to be set by the order of A and B
CC = B + A

# xarray Datasets are dict-like collections of DataArrays
# * dot indexing:  ds.salt
# * dict indexing: ds['salt'] - need to do this way for assignment, i.e. adding the variable,
# like when we added the z variables in extract_moor.py
