"""
Program to gather historical records for rivers.
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

# Load the as-run rivers from forcing files to fill in gaps in ec records
roms_fn = 'cas6_v3_2017.01.01_2020.12.31.p' # copied from LiveOcean_output/river
roms_qt = pd.read_pickle(Ldir['LOo'] / 'pre' / 'river' / 'Data_roms' / roms_fn)

# location for output
out_dir = Ldir['LOo'] / 'pre' / 'river' / 'Data_historical'
Lfun.make_dir(out_dir)

# set time range
dt0 = datetime(1980,1,1)
dt1 = datetime(2020,12,31)

testing = False
if testing:
    df = df.loc[['skagit', 'fraser']]
    dt0 = datetime(2017,1,1)
    dt1 = datetime(2020,12,31)

days = (dt0, dt1)
qt_dict = dict()
# get USGS river data
for rn in df.index:
    rs = df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.usgs):
        print('\n'+(' getting usgs ' + rn).center(60,'-'))
        rs, qt = rivf.get_usgs_data(rs, days)
        if rs['got_data']:
            qt_dict[rn] = qt
        else:
            print(rs)

# get ec data, a year at a time
for rn in df.index:
    rs = df.loc[rn].copy() # a series with info for this river
    Qt = pd.Series(dtype='float64') # initialize a Series to concatenate into
    if pd.notnull(rs.ec):# and rn in ['fraser']:
        for year in range(dt0.year, dt1.year + 1):
            this_days = (datetime(year,1,1), datetime(year,12,31))
            if year >= 2020:
                print((' ' + str(year) + ' getting ec ' + rn).center(60,'-'))
                rs, qt = rivf.get_ec_data(rs, this_days)
                
            # This section can be used as needed to fill gaps.
            # At this time 2021.04.09 there is no gap so we will omit.
            elif False: #year in [2017, 2018, 2019, 2020]:
                try:
                    qt = roms_qt.loc[this_days[0]:this_days[1], rn]
                    qt.index = qt.index + timedelta(days=0.5)
                    print((' ' + str(year) + ' getting ec from forcing ' + rn).center(60,'-'))
                    rs['got_data'] = True
                except Exception as e:
                    qt = ''
                    print(e)
                    rs['got_data'] = False
            
            elif year <= 2019:
                print((' ' + str(year) + ' getting ec historical ' + rn).center(60,'-'))
                rs, qt = rivf.get_ec_data_historical(rs, year)
                
            if rs.got_data:
                this_qt = qt.copy()
                Qt = pd.concat([Qt, this_qt])
            else:
                pass
                
        if not Qt.empty:
            qt_dict[rn] = Qt

# save output

for rn in qt_dict.keys():
    qt = qt_dict[rn]
    qt.to_pickle(out_dir / (rn + '.p'))

# plotting

if False:
    import matplotlib.pyplot as plt
    #plt.close('all')

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


