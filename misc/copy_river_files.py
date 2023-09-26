"""
This is a program to copy river forcing files from boiler onto apogee.

RESULT: I ran this on apogee on 2023.09.25 and it appeared to work.
Hence I should be able to simplify the river extraction code.
"""

from pathlib import Path
import pandas as pd
from lo_tools import Lfun
import shutil

in_dir0 = Path('/boildat/parker/LiveOcean_output/cas6_v3')

out_dir0 = Path('/dat1/parker/LO_output/forcing/cas6_v0')

date_range = pd.date_range(start='2016-12-15', end='2021-10-15')

f_list = ['f'+ item.strftime(Lfun.ds_fmt) for item in date_range]

for f in f_list:
    in_file = in_dir0 / f / 'riv2' / 'rivers.nc'
    out_dir = out_dir0 / f / 'riv0'
    Lfun.make_dir(out_dir)
    out_file = out_dir / 'rivers.nc'
    try:
        shutil.copy(in_file,out_file)
    except Exception as e:
        print(e)