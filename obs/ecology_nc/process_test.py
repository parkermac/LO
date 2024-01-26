"""
Code to explore the new NetCDF Ecology Archive.
"""

from lo_tools import Lfun
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

in_dir = Ldir['data'] / 'obs' / 'ecology_nc'
fn = in_dir / '1999to2023CTDnuts.nc'

ds = xr.open_dataset(fn, decode_times=False)

# converting the times
tu = ds.UTCDatetime.values
tu_dc = [item.decode() for item in tu]
ti = pd.to_datetime(tu_dc)

c = ds.CastGUID.values
c_dc = [item.decode() for item in c]

t = ds.Temp.values
p = ds.Pres.values
d = ds.DepthInterval.values

df = pd.DataFrame(index=ti, data={'c':c_dc,'t':t,'p':p})
