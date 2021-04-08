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

# location for output
out_dir = Ldir['LOo'] / 'pre' / 'river' / 'Data_historical'
Lfun.make_dir(out_dir)

# decide which group to get
get_usgs = True
get_ec = True
# and decide whether or not to save the data
save_data = True


# set time range
dt0 = datetime(1980,1,1)
dt1 = datetime(2020,12,31)

testing = True
if testing:
    df = df.loc[['skagit', 'fraser']]
    dt0 = datetime(2018,1,1)
    dt1 = datetime(2018,1,31)
    rs = df.loc['skagit'].copy()
    rs['flow_units'] = ''
    
days = (dt0, dt1)

if testing:
    rs, qt = rivf.get_usgs_data_sub(rs, days)

#
# qt_dict = dict()
#
# #%% get USGS river data
#
# if get_usgs:
#     for rn in df.index:
#         rs = df.loc[rn] # a series with info for this river
#         riv = river_class.River(rn, rs)
#         if pd.notnull(rs.usgs):
#             riv.get_usgs_data(days)
#             riv.print_info()
#             sys.stdout.flush()
#         if not riv.qt.empty:
#             qt_dict[rn] = riv.qt
#
# #%% get EC data, a year at a time
# #roms_fn = 'cas4_v2_2017.01.01_2018.12.31.p'
# roms_fn = 'cas6_v3_2017.01.01_2020.12.31.p'
# roms_qt = pd.read_pickle(Ldir['LOo'] + 'river/' + roms_fn)
# if get_ec:
#     for rn in df.index:
#         rs = df.loc[rn] # a series with info for this river
#         Qt = pd.Series() # initialize a Series to concatenate into
#         if pd.notnull(rs.ec):# and rn in ['fraser']:
#             for year in range(dt0.year, dt1.year + 1):
#                 print('year = ' + str(year))
#                 this_days = (datetime(year,1,1), datetime(year,12,31))
#                 riv = river_class.River(rn, rs)
#                 if year >= 2020:
#                     print(' - getting current EC data')
#                     riv.get_ec_data(this_days)
#                     riv.print_info()
#                     sys.stdout.flush()
#                     this_qt = riv.qt
#                 elif year in [2017, 2018, 2019, 2020]:
#                     this_qt = roms_qt.loc[this_days[0]:this_days[1], rn]
#                     this_qt.index = this_qt.index + timedelta(days=0.5)
#                     print(' - getting historical EC data from ' + roms_fn)
#                 elif year <= 2016:
#                     print(' - getting historical EC data')
#                     riv.get_ec_data_historical(year)
#                     riv.print_info()
#                     sys.stdout.flush()
#                     this_qt = riv.qt
#                 Qt = pd.concat([Qt, this_qt])
#             if not Qt.empty:
#                 qt_dict[rn] = Qt
#
# #%% save output
#
# if save_data:
#     for rn in qt_dict.keys():
#         qt = qt_dict[rn]
#         qt.to_pickle(out_dir + rn + '.p')
# else:
#     print('Not saving any data.')
#
# #%% plotting
#
# if True:
#
#     plt.close('all')
#
#     NP = len(qt_dict)
#     NR, NC = zfun.get_rc(NP)
#     fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
#     ii = 0
#     for rn in qt_dict.keys():
#         ir, ic = zfun.get_irc(ii, NC)
#         ax = axes[ir, ic]
#         this_ser = qt_dict[rn]
#         this_ser.plot(ax=ax, style='-k')
#         ax.set_xlim(dt0, dt1)
#         ax.text(.05, .9, rn, transform=ax.transAxes)
#         ii += 1
#
#     plt.show()


