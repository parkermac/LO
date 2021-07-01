"""
Plot bulk fluxes as a time series, including all chemical tracers
if present.

"""
from pathlib import Path
import sys

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import matplotlib.pyplot as plt
import numpy as np
import pickle

import Lfun
import zfun
import plotting_functions as pfun
import tef_fun
import flux_fun

Ldir = Lfun.Lstart()

in_dir00 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir00)
in_dir0 = in_dir00 / gtagex / 'tef'
ext_name = Lfun.choose_item(in_dir0, tag='bulk', exclude_tag='bulk_plots')
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.p')]
    
out_dir = in_dir0 / ext_name.replace('bulk', 'bulk_plots')
Lfun.make_dir(out_dir, clean=True)

# =================================
# get the DataFrame of all sections
gridname=gtagex.split('_')[0]
sect_df = tef_fun.get_sect_df(gridname)

testing = False

if testing:
    from importlib import reload
    reload(flux_fun)

if testing:
    sect_list = ['jdf1']
else:
    sect_list = list(sect_df.index)

# PLOTTING
fs = 14
pfun.start_plot(fs=fs, figsize=(17,10))

for sect_name in sect_list:

    # ---------------------------------------------------------

    tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(in_dir, sect_name, gridname)
    cols = list(tef_df.columns)
    
    # Check if we have chemical variables, or just salt
    vn_list = tef_fun.vn_list
    do_chem = True
    for vn in vn_list:
        if vn+'_in' not in cols:
            do_chem = False
    if do_chem:
        NR = 3; NC = 2 # rows and columns for subplots
    else:
        NR = 2; NC = 1
    
    # some information about direction
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        if in_sign == 1:
            dir_str = 'Eastward'
        elif in_sign == -1:
            dir_str = 'Westward'
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        if in_sign == 1:
            dir_str = 'Northward'
        elif in_sign == -1:
            dir_str = 'Southward'
                
    fig = plt.figure()

    # transport axis limits
    qlim_p = np.around(1.2*tef_df['Qin'].max(), 0)
    qlim_m = np.around(-1.2*tef_df['Qout'].min(), 0)
    qlim = np.max([qlim_p, qlim_m])
    
    # transport
    ax = fig.add_subplot(NR, NC, 1)
    tef_df[['Qin','Qout']].plot(ax=ax, legend=False, color=['r','b'], grid=True)
    ax.set_title(gtagex + ' : ' + sect_name + ' : Positive is ' + dir_str)
    ax.set_ylim(-qlim, qlim)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Transport $[m^{3}s^{-1}]$')
    ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', color='r',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', color='b',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
    # all other variables
    if do_chem:
        vn_list = ['salt', 'temp', 'NO3', 'oxygen', 'TIC']
        ylab_list = [r'Salinity', r'Pot. Temp. $[^{o}C]$', r'Nitrate $[\mu mol\ L^{-1}]$',
            r'Oxygen $[\mu mol\ L^{-1}]$', 'TIC (solid), Alk. (dashed)']
        ii_list = [3, 5, 2, 4, 6]
        xax_list = [False, True, False, False, True]
        abc_list = ['c', 'e', 'b', 'd', 'f']
        ii_dict = dict(zip(vn_list, ii_list))
        xax_dict = dict(zip(vn_list, xax_list))
        abc_dict = dict(zip(vn_list, abc_list))
        ylab_dict = dict(zip(vn_list, ylab_list))
    else:
        vn_list = ['salt']
        ylab_list = [r'Salinity']
        ii_list = [2]
        xax_list = [True]
        abc_list = ['b']
        ii_dict = dict(zip(vn_list, ii_list))
        xax_dict = dict(zip(vn_list, xax_list))
        abc_dict = dict(zip(vn_list, abc_list))
        ylab_dict = dict(zip(vn_list, ylab_list))
        
    for vn in vn_list:
        ax = fig.add_subplot(NR, NC, ii_dict[vn])
        tef_df[[vn+'_in',vn+'_out']].plot(ax=ax, legend=False, color=['r','b'], grid=True)
        if vn == 'TIC':
            tef_df[['alkalinity_in','alkalinity_out']].plot(ax=ax, legend=False,
                color=['r','b'], linestyle='--', grid=True)
        ax.set_ylabel(ylab_dict[vn])
        if xax_dict[vn]:
            pass
        else:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        ax.text(.03, .95, '(' + abc_dict[vn] + ')', va='top', weight='bold',
            transform=ax.transAxes, size=1.2*fs,
            bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
                
    fig.tight_layout()
    
    if testing:
        plt.show()
    else:
        plt.savefig(out_dir / (sect_name + '.png'))
        plt.close()

pfun.end_plot()
