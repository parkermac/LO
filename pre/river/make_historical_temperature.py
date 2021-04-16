"""
Program to gather historical records for river temperature.
"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import zfun
Ldir = Lfun.Lstart()

import river_functions as rivf
from importlib import reload
reload(rivf)

from datetime import datetime, timedelta
import pandas as pd
import numpy as np

# Load a dataframe with info for rivers to get
ri_fn = Ldir['LOo'] / 'pre' / 'river' / 'pnw_all_2021_04' / 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

# location for output
out_dir = Ldir['LOo'] / 'pre' / 'river' / 'Data_historical'
Lfun.make_dir(out_dir)

testing = False
if testing:
    df = df.loc[['clowhom']]

qt_dict = dict()

if True:
    # get USGS river data
    dt0 = datetime(1980,1,1)
    dt1 = datetime(2020,12,31)
    days = (dt0, dt1)
    for rn in df.index:
        rs = df.loc[rn].copy() # a series with info for this river
        if pd.notnull(rs.usgs):
            print((' getting usgs ' + rn).center(60,'-'))
            rs, qt = rivf.get_usgs_data(rs, days, temperature=True)
            if rs.got_data:
                qt_dict[rn] = qt
            else:
                print(rs)

if True:
    # get ec data
    dt0 = datetime(2020,1,1)
    dt1 = datetime(2020,12,31)
    days = (dt0, dt1)
    for rn in df.index:
        rs = df.loc[rn].copy() # a series with info for this river
        if pd.notnull(rs.ec):
            print((' getting ec ' + rn).center(60,'-'))
            rs, qt = rivf.get_ec_data(rs, days, temperature=True)
            if rs.got_data:
                qt[qt>100] = np.nan
                qt_dict[rn] = qt
            else:
                print(rs)

# save output

for rn in qt_dict.keys():
    qt = qt_dict[rn]
    if not testing:
        qt.to_pickle(out_dir / (rn + '_temperature.p'))

# plotting

if testing:
    import matplotlib.pyplot as plt
    plt.close('all')

    NP = len(qt_dict)
    NR, NC = zfun.get_rc(NP)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
    ii = 0
    for rn in qt_dict.keys():
        ir, ic = zfun.get_irc(ii, NC)
        ax = axes[ir, ic]
        this_ser = qt_dict[rn]
        this_ser.plot(ax=ax, style='-k')
        ax.set_xlim(dt0, dt1)
        ax.text(.05, .9, rn, transform=ax.transAxes)
        ii += 1

    plt.show()


