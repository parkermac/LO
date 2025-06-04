"""
Code to process the Willapa Bay mooring data.

"""

import pandas as pd
import xarray as xr
import numpy as np
import gsw
from time import time
from datetime import datetime, timedelta

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

# input location
source = 'willapa'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle

year_list = [2017]

if testing:
    from lo_tools import plotting_functions as pfun
    import matplotlib.pyplot as plt
    plt.close('all')
    pfun.start_plot()

for year in year_list:
    year_str = str(year)

    in_dir = Ldir['data'] / 'obs' / source / otype / year_str

    # output location
    out_dir = Ldir['LOo'] / 'obs' / source / otype / year_str
    if not testing:
        Lfun.make_dir(out_dir, clean=True)

    fn_dict = {
        'BayCenter':'Bay Center 2017share.xlsx', # is this correct?
        'Nahcotta':'Nahcotta 2017share.xlsx',
    }

    sta_dict = {
        'Tokeland': (-123.96797728985608, 46.707252252710354),
        'Nahcotta': (-124.03074795328776, 46.500242867945865),
        'BayCenter': (-123.95239473341415, 46.629030151420984),
    }

    if False:
        sta_list = ['Tokeland']
    else:
        sta_list = list(fn_dict.keys())

    for sta in sta_list:
        in_fn = in_dir / fn_dict[sta]
        out_fn = out_dir / (sta + '.p') # a pickled DataFrame
    
        a = pd.read_excel(in_fn)

        # "ti" local time, but when I tried to carefully convert to UTC
        # I ran into this error:
        # NonExistentTimeError: 2017-03-12 02:00:00
        # which is because the shift to PDT skips this hour.
        # I will proceed by just adding 8 hours to make it
        # "UTC without the Daylight Savings shift"
        # 
        # For future reference, this is the code I tried to use
        # import pytz
        # tz = pytz.timezone('America/Los_Angeles')
        # # Localize the DatetimeIndex to the local timezone
        # til = ti.tz_localize(tz)
        # # Convert to UTC
        # tiu = til.tz_convert('UTC')

        if sta == 'BayCenter':
            ti = pd.DatetimeIndex(a['Date'])
            SP = a['S Corr'].to_numpy()
            IT = a['Temp C'].to_numpy()
            PH = a['pH Corr'].to_numpy()
            do_ph = True
        elif sta == 'Nahcotta':
            ti = pd.DatetimeIndex(a['Date Time'])
            SP = a['Sal Corr'].to_numpy()
            IT = a['T ©'].to_numpy()
            do_ph = False
        lon, lat = sta_dict[sta]
        P = 0 # pressure (the mooring is near the surface)
        # - do the conversions
        SA = gsw.SA_from_SP(SP, P, lon, lat)
        CT = gsw.CT_from_t(SA, IT, P)

        # pack in a DataFrame
        tiu = ti + timedelta(hours=8)
        if do_ph:
            df = pd.DataFrame(index=tiu, columns=['SA', 'CT', 'PH'])
            df['SA'] = SA
            df['CT'] = CT
            df['PH'] = PH
        else:
            df = pd.DataFrame(index=tiu, columns=['SA', 'CT'])
            df['SA'] = SA
            df['CT'] = CT

        if testing:
            # for some reason this produces 4 figures, 2 of them blank
            fig = plt.figure()
            df.plot(fig=fig, subplots=True, title=sta)
        else:
            df.to_pickle(out_fn)

