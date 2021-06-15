"""
Plot the exchange flow in a dynamical context, using the time mean
of all sections

"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from time import time

import Lfun
import zfun
import plotting_functions as pfun
import tef_fun
import flux_fun

Ldir = Lfun.Lstart()

gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)

ii = 0
for gtagex in ['cas6_v3_lo8b', 'cas6_v3t075_lo8']:
    in_dir = Path('/Users/pm8/Documents/LO_output/extract/'+gtagex+'/tef/bulk_2018.01.01_2018.12.31')
    pro_dir = Path(str(in_dir).replace('bulk','processed'))
    tt0 = time()
    df = pd.DataFrame(index=sect_list)
    for sect_name in sect_list:
        # get results of process_sections
        pro = pickle.load(open(pro_dir / (sect_name + '.p'), 'rb'))
        qnet = pro['qnet'].copy() # hourly time series
        qabs = np.abs(qnet)
        pad = 34
        qtide = (np.pi/2) * zfun.filt_godin(qabs)[pad:-pad+1:24]
        qprism = qtide / np.pi
        # get two-layer time series
        gridname = 'cas6'
        tef_df, in_sign, dir_str = flux_fun.get_two_layer(in_dir, sect_name, gridname)
        # make derived variables
        df.loc[sect_name, 'Qe'] = ((tef_df['Qin'] - tef_df['Qout']).mean()/2)/1000
        df.loc[sect_name, 'DS'] = (tef_df['salt_in'] - tef_df['salt_out']).mean()
        df.loc[sect_name, 'Sbar'] = (tef_df['salt_in'] + tef_df['salt_out']).mean()/2
        df.loc[sect_name, 'Qprism'] = np.nanmean(qprism)/1000
    if ii == 0:
        df0 = df.copy()
    elif ii == 1:
        df1 = df.copy()
    print('Total time to fill DataFrame = %0.1f sec' % (time()-tt0))
    ii += 1
    
dfx = pd.DataFrame(index=sect_list)
dfx['r0'] = df0['Qe']/df0['Qprism']
dfx['r1'] = df1['Qe']/df1['Qprism']
ax = dfx.plot(x='r0', y='r1', style='*g')
ax.plot([0,1],[0,1],'-k')

# # PLOTTING
# ax = df0.plot(x='Qprism', y='Qe', style='*b')
# df1.plot(x='Qprism', y='Qe', style='*r', ax=ax)
