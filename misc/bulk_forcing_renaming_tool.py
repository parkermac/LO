"""
This custom code will rename a series of forcing folders over
a date range. For example if you wanted to rename ocnG01 -> ocnG02
in all folders in LO_output/forcing/cas7 over some date range.
"""

from lo_tools import Lfun
import pandas as pd
import os
from datetime import datetime

Ldir = Lfun.Lstart()
gridname = 'cas7'
old_name = 'ocn03'
new_name = 'ocn03_test'
dt0 = datetime(2024,11,22)
dt1 = datetime(2024,11,23)

dti = pd.date_range(dt0,dt1)

for dt in dti:
    print(dt.strftime(Ldir['ds_fmt']))
    fstring = 'f' + dt.strftime(Ldir['ds_fmt'])
    old_dir = Ldir['LOo'] / 'forcing' / gridname / fstring / old_name
    new_dir = Ldir['LOo'] / 'forcing' / gridname / fstring / new_name
    os.rename(old_dir,new_dir)
