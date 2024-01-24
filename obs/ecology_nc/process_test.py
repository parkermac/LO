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
t = ds.UTCDatetime.values
tt = [item.decode() for item in t]
ti = pd.to_datetime(tt)