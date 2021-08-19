"""
Plot the mean of many TEF extractions on all the channels.

Clean version for publication.

"""
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun
import flux_fun

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df(gridname)

# select input directory
in_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
in_fn = in_dir / 'two_layer_mean_2018.08.01_2018.12.31.p'

def plotit(ax, sect_df, sect_list, lcol, qsign, alpha=.5):
    counter = 0
    for sn in sect_list:
        # some information about direction
        x0, x1, y0, y1 = sect_df.loc[sn,:]
        ax.plot([x0,x1], [y0,y1], '-', color=lcol, lw=3, alpha=alpha)
        xx = (x0+x1)/2
        yy = (y0+y1)/2
        # add a mark showing which direction the deep flow is going,
        # noting that the transports are now all positive east or north.
        clat = np.cos(np.pi*yy/180)
        if (x0==x1) and (y0!=y1):
            sdir = 'NS'
            dd = qsign[counter] * 0.05 / clat
            ww = dd/4
            ax.fill([xx, xx+dd, xx], [yy-ww, yy, yy+ww], color=lcol, alpha=alpha)
        elif (x0!=x1) and (y0==y1):
            sdir = 'EW'
            dd = qsign[counter] * 0.05
            ww = dd/(4*clat)
            ax.fill([xx-ww, xx, xx+ww], [yy, yy+dd, yy], color=lcol, alpha=alpha)
        counter += 1
    
def inax(x,y,aa):
    # returns True if x,y, is inside of the box defined by aa
    is_in = x>aa[0] and x<aa[1] and y>aa[2] and y<aa[3]
    return is_in

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
channel_list.reverse() # makes overlaying colors look better

distmax = 420

aa1 = [-5, distmax, 0 ,200]
aa2 = [140, 300, -5, 50]
aa3 = [-125.4, -122, 46.8, 50.4] # Salish Sea

# load the two_layer_mean DataFrame
df = pd.read_pickle(in_fn)

fig = plt.figure(figsize=(14,10))

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(122)

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
    print('== ' + ch_str + ' ==')
    sect_list = channel_dict[ch_str]

    q_s = df.loc[sect_list,'Qin'].to_numpy(dtype='float')/1e3
    q_f = df.loc[sect_list,'Qout'].to_numpy(dtype='float')/1e3
    s_s = df.loc[sect_list,'salt_in'].to_numpy(dtype='float')
    s_f = df.loc[sect_list,'salt_out'].to_numpy(dtype='float')
    # get the sign of q_s and plot the section locations with "inflow" direction
    # (defined as the direction of transport of the saltier water of the pair)
    qsign = df.loc[sect_list,'in_sign'].to_numpy(dtype='float')#np.sign(q_s)
    dist = dist_dict[ch_str]

    lcol = lcol_dict[ch_str]

    if ch_str != 'Juan de Fuca to Strait of Georgia':
        ax1.plot(dist[1:],np.abs(q_s[1:]),'-', color=lcol, lw=lw)
        ax1.plot(dist[1:],np.abs(q_f[1:]),'-', color=lcol, lw=lw)
        ax1.plot(dist[:2],np.abs(q_s[:2]),'--', color='k', lw=1)
        ax1.plot(dist[:2],np.abs(q_f[:2]),'--', color='k', lw=1)
    else:
        ax1.plot(dist,np.abs(q_s),'-', color=lcol, lw=lw)
        ax1.plot(dist,np.abs(q_f),'-', color=lcol, lw=lw)

    counter = 0
    for sn in sect_list:
        xx = dist[counter]
        yy = np.abs(q_s[counter])
        if inax(xx,yy,aa1) and (ch_str == 'Juan de Fuca to Strait of Georgia'):
            dy = aa1[3]-aa1[2]
            ax1.text(xx,yy + dy/10, sn, fontsize=.75*fs, color='k', va='center',ha='center',
                rotation=45, style='italic', weight='bold')
        counter += 1
    
    lcol = lcol_dict[ch_str]

    if ch_str != 'Juan de Fuca to Strait of Georgia':
        ax2.plot(dist[1:],np.abs(q_s[1:]),'-', color=lcol, lw=lw)
        ax2.plot(dist[1:],np.abs(q_f[1:]),'-', color=lcol, lw=lw)
        ax2.plot(dist[:2],np.abs(q_s[:2]),'--', color='k', lw=1)
        ax2.plot(dist[:2],np.abs(q_f[:2]),'--', color='k', lw=1)
    else:
        ax2.plot(dist,np.abs(q_s),'-', color=lcol, lw=lw)
        ax2.plot(dist,np.abs(q_f),'-', color=lcol, lw=lw)

    counter = 0
    for sn in sect_list:
        xx = dist[counter]
        yy = np.abs(q_s[counter])
        if inax(xx,yy,aa2):
            dy = aa2[3]-aa2[2]
            if ch_str == 'Hood Canal' and sn != 'ai3':
                ax2.text(xx,yy - dy/20, sn, fontsize=.75*fs, color='k', va='center',ha='center',
                    rotation=-45, style='italic', weight='bold')
            else:
                ax2.text(xx,yy + dy/20, sn, fontsize=.75*fs, color='k', va='center',ha='center',
                    rotation=45, style='italic', weight='bold')
        
        counter += 1
    ax2.set_xlabel('Distance from Mouth [km]')

    # MAP
    plotit(ax3, sect_df, sect_list, lcol, qsign)
    

ax1.axis(aa1)
ax1.grid(True)
ax1.set_ylabel('Transport $[1000\ m^{3}s^{-1}]$')
#ax1.set_title(season.title() + ' ' + year_str)

ax1.fill([aa2[0], aa2[1], aa2[1], aa2[0], aa2[0]], [aa2[2], aa2[2], aa2[3], aa2[3], aa2[2]],
    'orange', alpha=.3)

ax2.axis(aa2)
ax2.grid(True)
ax2.set_ylabel('Transport $[1000\ m^{3}s^{-1}]$')

pfun.add_coast(ax3, color='gray')
pfun.dar(ax3)
ax3.axis(aa3)
ax3.set_xticks([])
ax3.set_yticks([])
# just do these things once
for sn in sect_df.index:
    x0, x1, y0, y1 = sect_df.loc[sn,:]
    ax3.text(x1,y1, sn, fontsize=.6*fs, color='k', va='center',ha='center',
        rotation=45, style='italic', weight='bold')

ax1.text(.97, .95, '(a)', ha='right', va='top', weight='bold', transform=ax1.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax2.text(.97, .95, '(b)', ha='right', va='top', weight='bold', transform=ax2.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax3.text(.97, .95, '(c)',
    ha='right', va='top', weight='bold', transform=ax3.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax3.text(.97, .9, 'Section Locations\n& Deep Inflow Directions',
    ha='right', va='top', transform=ax3.transAxes, style='italic',
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

fig.tight_layout()
plt.show()
#fig.savefig(outdir + 'all_sections_' + year_str + '_' + season + '.png')
    
plt.rcdefaults()
    