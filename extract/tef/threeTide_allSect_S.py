"""
Plot the two-layer mean salinity of many TEF extractions on all the channels.

Compares three different runs with varying tidal strength.

"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import Lfun
Ldir = Lfun.Lstart()

import plotting_functions as pfun

import tef_fun
import flux_fun

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df('cas6')

# PLOTTING
lw=3
fs=16
ms = 20
alpha = .5
qscl = 50
plt.rc('font', size=fs)
plt.close('all')

channel_list = flux_fun.channel_list
channel_dict = flux_fun.long_channel_dict

lcol_dict = flux_fun.c_dict

ch_idict = dict(zip(channel_list, [1,2,3,4])) # subplot number

# channel_list =
# ['Whidbey Basin',
#  'Hood Canal',
#  'Admiralty Inlet to South Sound',
#  'Juan de Fuca to Strait of Georgia']
aad = dict()
aad[1] = [190, 260, 23, 32]
aad[2] = [170, 270, 28.5, 32]
aad[3] = [140, 300, 30, 32.5]
aad[4] = [0, 400, 29, 33.5]

fig = plt.figure(figsize=(14,10))

axd = dict()
axd[1] = fig.add_subplot(221)
axd[2] = fig.add_subplot(222)
axd[3] = fig.add_subplot(223)
axd[4] = fig.add_subplot(224)

gtagex_list = ['cas6_v3t075_lo8', 'cas6_v3_lo8b', 'cas6_v3t110_lo8']
alpha_dict = dict(zip(gtagex_list, [.3, .6, 1]))

# loop over all three gtagex
for gtagex in gtagex_list:
    
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef'
    in_fn = in_dir / 'two_layer_mean_2018.08.01_2018.12.31.p'

    # load the two_layer_mean DataFrame
    df = pd.read_pickle(in_fn)

    # create all the distance vectors and save in a dict
    dist_dict = {}
    for ch_str in channel_list:
        sect_list = channel_dict[ch_str]
        x = df.loc[sect_list,'lon'].to_numpy(dtype='float')
        y = df.loc[sect_list,'lat'].to_numpy(dtype='float')
        dist_dict[ch_str] = flux_fun.make_dist(x,y)

    # adjust the distance vectors to join at the correct locations
    ind_ai = channel_dict['Juan de Fuca to Strait of Georgia'].index('jdf4')
    dist0_ai = dist_dict['Juan de Fuca to Strait of Georgia'][ind_ai]
    dist_dict['Admiralty Inlet to South Sound'] += dist0_ai
    #
    ind_hc = channel_dict['Admiralty Inlet to South Sound'].index('ai3')
    dist0_hc = dist_dict['Admiralty Inlet to South Sound'][ind_hc]
    dist_dict['Hood Canal'] += dist0_hc
    #
    ind_wb = channel_dict['Admiralty Inlet to South Sound'].index('ai4')
    dist0_wb = dist_dict['Admiralty Inlet to South Sound'][ind_wb]
    dist_dict['Whidbey Basin'] += dist0_wb

    for ch_str in channel_list:
        sect_list = channel_dict[ch_str]
        s_s = df.loc[sect_list,'salt_in'].to_numpy(dtype='float')
        s_f = df.loc[sect_list,'salt_out'].to_numpy(dtype='float')
        qsign = df.loc[sect_list,'in_sign'].to_numpy(dtype='float')
        dist = dist_dict[ch_str]

        lcol = lcol_dict[ch_str]
        
        axd[ch_idict[ch_str]].plot(dist[1:],(s_s[1:]),'-',
            color=lcol, lw=lw, alpha=alpha_dict[gtagex])
        axd[ch_idict[ch_str]].plot(dist[1:],(s_f[1:]),'--',
            color=lcol, lw=lw,alpha=alpha_dict[gtagex])
            
    
for ch_str in channel_list:
    lcol = lcol_dict[ch_str]
    ax = axd[ch_idict[ch_str]]
    ax.text(.05,.05,ch_str,c=lcol, transform=ax.transAxes, weight='bold')
    
for ii in [1,2,3,4]:
    axd[ii].axis(aad[ii])
    
fig.tight_layout()
plt.show()
#fig.savefig(outdir + 'all_sections_' + year_str + '_' + season + '.png')
    
plt.rcdefaults()
    