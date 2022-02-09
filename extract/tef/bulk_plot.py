"""
Plot bulk fluxes as a time series, including all chemical tracers
if present.

"""
import matplotlib.pyplot as plt
import numpy as np
import pickle

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
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
fs = 12
plt.close('all')
pfun.start_plot(fs=fs, figsize=(21,10))

for sect_name in sect_list:

    tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(in_dir, sect_name, gridname)
    cols = list(tef_df.columns)
    
    # Check if we have chemical variables, or just salt
    vn_list = tef_fun.vn_list
    do_chem = True
    for vn in vn_list:
        if vn+'_in' not in cols:
            do_chem = False
            
    if do_chem:
        vn_list = vn_list = ['Q', 'salt', 'temp', 'NO3', 'phytoplankton', 'detritus', 'oxygen', 'TIC', 'alkalinity']
        vn_list_long = vn_list + ['zooplankton', 'Ldetritus']
        NR = 3; NC = 3 # rows and columns for subplots
    else:
        vn_list = ['Q', 'salt']
        vn_list_long = vn_list
        NR = 2; NC = 1
        
    # filter more in time (requires that Q be the first item in vn_list)
    nhan = 1 # length of Hanning window in days (use 1 for no filtering)
    for vn in vn_list_long:
        if vn == 'Q':
            tef_df['Q_in'] = zfun.lowpass(tef_df['Qin'].to_numpy(), n=nhan)
            tef_df['Q_out'] = zfun.lowpass(tef_df['Qout'].to_numpy(), n=nhan)
        else:
            tef_df[vn+'_in'] = zfun.lowpass((tef_df[vn+'_in']*tef_df['Qin']).to_numpy(), n=nhan)/tef_df['Q_in']
            tef_df[vn+'_out'] = zfun.lowpass((tef_df[vn+'_out']*tef_df['Qout']).to_numpy(), n=nhan)/tef_df['Q_out']
            
    # adjust units
    tef_df['Q_in'] = tef_df['Q_in']/1000
    tef_df['Q_out'] = tef_df['Q_out']/1000
    
                    
    # labels and colors
    ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
                'salt': r'Salinity',
                'temp': r'Pot. Temp. $[^{o}C]$',
                'NO3': r'Nitrate $[\mu mol\ N\ L^{-1}]$',
                'phytoplankton': r'phytoplankton $[\mu mol\ N\ L^{-1}]$',
                'detritus': r'detritus $[\mu mol\ N\ L^{-1}]$',
                'oxygen': r'oxygen $[\mu mol\ O_{2}\ L^{-1}]$',
                'TIC': r'TIC $[\mu mol\ C\ L^{-1}]$',
                'alkalinity': r'Alkalinity $[\mu equiv\ L^{-1}]$'}
    import string
    abc = list(string.ascii_lowercase)
    incolor = 'r'
    outcolor = 'b'
    lw = 2
    
    fig, axes = plt.subplots(nrows=NR, ncols=NC, squeeze=False)
    
    ii = 0
    for vn in vn_list:
        ir, ic = zfun.get_irc(ii, NC)
        ax = axes[ir, ic]
        tef_df[[vn+'_in',vn+'_out']].plot(ax=ax, legend=False, color=[incolor, outcolor], grid=True, lw=lw)
        ax.set_ylabel(ylab_dict[vn])
        ax.text(.03, .95, '('+abc[ii]+')', va='top', weight='bold', transform=ax.transAxes,
            bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        if ir < NR-1:
            ax.set_xticklabels([])
        if vn == 'Q':
            ax.set_title(gtagex + ' : ' + sect_name + ' : Positive is ' + dir_str)
            ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', c=incolor,
                transform=ax.transAxes, size=1.2*fs,
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
            ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', c=outcolor,
                transform=ax.transAxes, size=1.2*fs,
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        if vn == 'phytoplankton':
            tef_df[['zooplankton_in','zooplankton_out']].plot(ax=ax, legend=False, color=[incolor, outcolor],
                grid=True, alpha=.5, lw=lw)
            ax.text(.97, .95, 'Transparent = zooplankton', ha='right', va='top', weight='bold', c='k',
                transform=ax.transAxes, size=fs, alpha=.5,
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        if vn == 'detritus':
            tef_df[['Ldetritus_in','Ldetritus_out']].plot(ax=ax, legend=False, color=[incolor, outcolor],
                grid=True, alpha=.5, lw=lw)
            ax.text(.97, .95, 'Transparent = Ldetritus', ha='right', va='top', weight='bold', c='k',
                transform=ax.transAxes, size=fs, alpha=.5,
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        ii += 1
                
    fig.tight_layout()
    
    if testing:
        plt.show()
    else:
        plt.savefig(out_dir / (sect_name + '.png'))
        plt.close()

pfun.end_plot()
