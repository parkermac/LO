"""
Code to create the mean over a specified time period of the two-layer
TEF properties, and some derived quantities, for several runs.

It takes 9 sec or so for each gtagex.

The names and time periods are hard-coded, so you have to edit it for
different purposes.

"""
from pathlib import Path
import sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from time import time

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun
import flux_fun

Ldir = Lfun.Lstart()

gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)

ds0 = '2018.01.01'
ds1 = '2018.12.31'
dt00 = datetime(2018,1,1)
ds00 = dt00.strftime(Lfun.ds_fmt)
ii = 0
for gtagex in ['cas6_v3_lo8b', 'cas6_v3t075_lo8', 'cas6_v3t110_lo8']:
    in_dir = Path('/Users/pm8/Documents/LO_output/extract/'+gtagex+'/tef/bulk_' + ds0 + '_' + ds1)
    tt0 = time()
    df = pd.DataFrame(index=sect_list)
    for sect_name in sect_list:
        x0, x1, y0, y1 = sect_df.loc[sect_name,:]
        lon = (x0 + x1)/2
        lat = (y0 + y1)/2
        # get two-layer versions
        tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(in_dir, sect_name, gridname, dt00=dt00)
        # and save their flux-weighted means
        df.loc[sect_name, 'Qin'] = tef_df['Qin'].mean()
        df.loc[sect_name, 'Qout'] = tef_df['Qout'].mean()
        for vn in tef_fun.vn_list:
            if (vn+'_in') in tef_df.columns:
                df.loc[sect_name, vn+'_in'] = (tef_df['Qin']*tef_df[vn+'_in']).mean()/tef_df['Qin'].mean()
                df.loc[sect_name, vn+'_out'] = (tef_df['Qout']*tef_df[vn+'_out']).mean()/tef_df['Qout'].mean()
        # make derived variables
        Socn = 34
        df.loc[sect_name, 'Qfw'] = (tef_df['Qin']*(Socn-tef_df['salt_in']) + tef_df['Qout']*(Socn-tef_df['salt_out'])).mean()/Socn
        df.loc[sect_name, 'Qe'] = ((df.loc[sect_name, 'Qin'] - df.loc[sect_name, 'Qout'])/2)/1000
        df.loc[sect_name, 'DS'] = df.loc[sect_name, 'salt_in'] - df.loc[sect_name, 'salt_out']
        df.loc[sect_name, 'Sbar'] = (df.loc[sect_name, 'salt_in'] + df.loc[sect_name, 'salt_out'])/2
        df.loc[sect_name, 'Qprism'] =( tef_df['qabs'].mean()/2)/1000
        df.loc[sect_name, 'SSH'] = tef_df['ssh'].mean()
        df.loc[sect_name, 'Fnet'] = tef_df['fnet'].mean()
        df.loc[sect_name, 'Qnet'] = tef_df['qnet'].mean()
        df.loc[sect_name, 'in_sign'] = in_sign
        df.loc[sect_name, 'dir_str'] = dir_str
        df.loc[sect_name, 'sdir'] = sdir
        df.loc[sect_name, 'lon'] = lon
        df.loc[sect_name, 'lat'] = lat
    print('Total time to fill DataFrame = %0.1f sec' % (time()-tt0))
    ii += 1
    df.to_pickle(in_dir.parent / ('two_layer_mean_' + ds00 + '_' + ds1 + '.p'))


