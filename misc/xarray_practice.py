"""
Code to practice and test the use of xarray in the LO system.
"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
import Lfun
import zrfun

import numpy as np
import pandas as pd
import xarray as xr

# open a ROMS history file
Ldir = Lfun.Lstart()
fn = (Ldir['parent'] / 'LiveOcean_roms' / 'output' /
    'cas6_v3_lo8b' / 'f2019.07.04' / 'ocean_his_0001.nc')
G,S,T = zrfun.get_basic_info(fn)

xs = xr.open_dataset(fn)
"""
An xarray Dataset with attibutes like:
- attrs
- data_vars (like salt)
- dims (like xi_rho)
- coords (lke lon_rho)
- values (typically an ndarray)

xr.open_dataset() = lazy loading
xr.load_dataset() = load everything
"""

s = xs.salt
# an xarray DataArray of shape (1, 30, 1302, 663)

# NOTE: hereafter we refer to a DataArray as simply an array

# make a DataArray from scratch (could also just pass a pandas Series for a 1D array)
a = xr.DataArray(np.arange(12).reshape((3,4)), dims=('x','y'), coords={'x':[1,2,3], 'y':list('abcd')})
# ways to work with it (these generally return another DataArray)
a = a * 2 # math
am = a.mean(dim='y')

# There are FOUR ways to index the array.
# - these two ways CAN change values in the array:
a[0, :] = 77 # positional and by integer label, like numpy
a.loc[2,'b'] = 88 # loc or "location": positional and coordinate label, like pandas
# - whereas these two ways can only return a new array, not change values in the array
aa = a.isel(x=0) # isel or "integer select":  by dimension name and integer label
bb = a.sel(y='b') # sel or "select": by dimension name and coordinate label

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
