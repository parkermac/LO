"""
Plot bulk fluxes as a time series.

To test on mac:
run bulk_plot -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import flux_fun
from importlib import reload
reload(flux_fun)

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('bulk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.p')]
if Ldir['testing']:
    sect_list = ['mb8.p']


# PLOTTING
fs = 12
plt.close('all')
pfun.start_plot(fs=fs, figsize=(21,10))

for sect_name in sect_list:
    
    bulk = pickle.load(open(in_dir / sect_name, 'rb'))

    tef_df = flux_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
                    
    # labels and colors
    ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
                'salt': r'Salinity $[g\ kg^{-1}]$'}
    p_color = 'r'
    m_color = 'b'
    lw = 2
        
    fig = plt.figure()
    
    ax1 = plt.subplot2grid((2,4), (0,0), colspan=3) # Qin, Qout
    ax2 = plt.subplot2grid((2,4), (1,0), colspan=3) # Sin, Sout
    ax3 = plt.subplot2grid((1,4), (0,3)) # map
    
    tef_df[['Q_p','Q_m']].plot(ax=ax1, legend=False, color=[p_color, m_color], grid=True, lw=lw)
    ax1.set_ylabel(ylab_dict['Q'])
    
    tef_df[['salt_p','salt_m']].plot(ax=ax2, legend=False, color=[p_color, m_color], grid=True, lw=lw)
    ax2.set_ylabel(ylab_dict['salt'])
        
                
    fig.tight_layout()
    
    if Ldir['testing']:
        plt.show()
    else:
        plt.savefig(out_dir / (sect_name + '.png'))
        plt.close()

pfun.end_plot()
