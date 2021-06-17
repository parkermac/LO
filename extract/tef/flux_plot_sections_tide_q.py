"""
Plot the mean of tidal energy flux and volume transport at all
TEF sections.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)


plt.close('all')

for year in [2017, 2018, 2019]:
    year_str = str(year)

    # input
    run_name = Ldir['gtagex']+'_'+year_str+'.01.01_'+year_str+'.12.31'
    indir00 = Ldir['LOo'] + 'tef/'
    indir0 = indir00 + run_name + '/'
    indir = indir0 + 'bulk/'

    # output
    out_fn = indir00 + 'misc_figs_' + Ldir['gridname'] + '/Tide_and_Qnet_' + year_str + '.png'

    # colors
    clist = flux_fun.clist

    # angles along the thalweg that are nominally "landward"
    # meant to be consistent with landward sign in sect_df
    th_dict = {
        'jdf1': -15,
        'jdf2': -15,
        'jdf3': -5,
        'jdf4': 0,
        'sog1': 135,
        'sog2': 160,
        'sog3': 150,
        'sog4': 100,
        'sog5': 90,
        'sji1': 90,
        'sji2': 100,
        'dp': 0,
        'ai1': -45,
        'ai2': -60,
        'ai3': -60,
        'ai4': -30,
        'wb1': 80,
        'wb2': 135,
        'wb3': 90,
        'wb4': 30,
        'hc1': -80,
        'hc2': -120,
        'hc3': -150,
        'hc4': -135,
        'hc5': -130,
        'hc6': -90,
        'hc7': 20,
        'hc8': 30,
        'mb1': -90,
        'mb2': -100,
        'mb3': -75,
        'mb4': -80,
        'mb5': -135,
        'tn1': -80,
        'tn2': -90,
        'tn3': -100,
        'ss1': -135,
        'ss2': 135,
        'ss3': -135,
    }

    # get the DataFrame of all sections
    sect_df = tef_fun.get_sect_df()
    sect_list = list(sect_df.index)

    # omit some sections for clarity
    sect_list.remove('tn1')
    sect_list.remove('tn3')
    
    df = pd.DataFrame(index=sect_list)
    for sn in sect_list:
        bulk = pickle.load(open(indir + sn + '.p', 'rb'))
        sx0, sx1, sy0, sy1, landward = sect_df.loc[sn,:]
        sx = (sx0+sx1)/2; sy = (sy0+sy1)/2
        if (sx0==sx1) and (sy0!=sy1):
            sdir = 'NS'
        elif (sx0!=sx1) and (sy0==sy1):
            sdir = 'EW'
        df.loc[sn,'lon'] = sx
        df.loc[sn,'lat'] = sy
        df.loc[sn,'sdir'] = sdir
        df.loc[sn,'landward'] = landward
        df.loc[sn,'F'] = bulk['fnet_lp'].mean()/1e6 # MW
        df.loc[sn,'Q'] = bulk['qnet_lp'].mean() # 1000 m3/s
    
    # PLOTTING
    fig = plt.figure(figsize=(10,13))
    fs = 18
    cf = 'purple'
    cq = 'g'

    # axis limits
    x0 = -125.5; x1 = -122; y0 = 47; y1 = 50.5 # Salish Sea
    x00 = -123.3; x11 = -122.2; y00 = 47; y11 = 48.5 # Puget Sound
    aaS = [x0, x1, y0, y1]
    aaP = [x00, x11, y00, y11]

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    for sn in df.index:
        sx = df.loc[sn,'lon']
        sy = df.loc[sn,'lat']
        sdir = df.loc[sn,'sdir']
        landward = df.loc[sn,'landward']
        # landward is the sign to multiply by to know what the direction
        # of a positive flux is.  So if sdir = 'NS' and landward = 1,
        # then positive flux is to the East
        F = df.loc[sn,'F']
        Q = df.loc[sn,'Q']
        logF = np.log10(np.abs(F) + 1)
        logQ = np.log10(np.abs(Q) + 1)
        sgnF = np.sign(landward*F)
        sgnQ = np.sign(landward*Q)
        scl = 40 # a vector of length 1 will be 1/scl of the length of the y-axis
        if False:
            if sdir == 'EW':
                vf0 = 0; vf1 = sgnF
                vq0 = 0; vq1 = sgnQ
            elif sdir == 'NS':
                vf0 = sgnF; vf1 = 0
                vq0 = sgnQ; vq1 = 0
        else:
            sth = np.sin(np.pi*th_dict[sn]/180)
            cth = np.cos(np.pi*th_dict[sn]/180)
            if landward == 1:
                vf0 = cth*sgnF*logF; vf1 = sth*sgnF*logF
                vq0 = cth*sgnQ*logQ; vq1 = sth*sgnQ*logQ
            elif landward == -1:
                vf0 = -cth*sgnF*logF; vf1 = -sth*sgnF*logF
                vq0 = -cth*sgnQ*logQ; vq1 = -sth*sgnQ*logQ
                
        aq = .5
        ffs = .7*fs
        ax1.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
            linewidths=1, color=cf, alpha=aq)
        if sx < x00 or sy > y11:
            ax1.text(sx, sy+.04, str(int(np.abs(F))), ha='center', va='center', size=ffs,
                weight='bold', color=cf, alpha=1)

        if sx > x00 and sy < y11:
            ax2.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
                linewidths=1, color=cf, alpha=aq)
            ax2.text(sx, sy+.02, str(int(np.abs(F))), ha='center', va='center', size=ffs,
            weight='bold', color=cf, alpha=1)

        ax3.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
            linewidths=1, color=cq, alpha=aq)
        if sx < x00 or sy > y11:
            ax3.text(sx, sy+.04, str(int(np.abs(Q))), ha='center', va='center', size=ffs,
                weight='bold', color=cq, alpha=1)

        if sx > x00 and sy < y11:
            ax4.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
                linewidths=1, color=cq, alpha=aq)
            ax4.text(sx, sy+.02, str(int(np.abs(Q))), ha='center', va='center', size=ffs,
            weight='bold', color=cq, alpha=1)

    pfun.add_coast(ax1, color='gray')
    ax1.axis(aaS)
    pfun.dar(ax1)
    ax1.set_xticklabels([])
    ax1.tick_params(labelsize=.8*fs)
    ax1.text(.95,.9,'(a)', size=fs, transform=ax1.transAxes, ha='right', weight='bold')
    ax1.set_ylabel('Latitude', size=.8*fs)
    ax1.set_xticks([-125, -124, -123, -122])
    ax1.text(.05,.05,'Tidal energy Flux $[MW]$', size=.8*fs, transform=ax1.transAxes,
        weight='bold', color=cf)

    pfun.add_coast(ax2, color='gray')
    ax2.axis(aaP)
    pfun.dar(ax2)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.tick_params(labelsize=.8*fs)
    ax2.text(.95,.9,'(b)', size=fs, transform=ax2.transAxes, ha='right', weight='bold')
    ax2.set_xticks([-123, -122.5])

    pfun.add_coast(ax3, color='gray')
    ax3.axis(aaS)
    pfun.dar(ax3)
    ax3.tick_params(labelsize=.8*fs)
    ax3.text(.95,.9,'(c)', size=fs, transform=ax3.transAxes, ha='right', weight='bold')
    ax3.set_xlabel('Longitude', size=.8*fs)
    ax3.set_ylabel('Latitude', size=.8*fs)
    ax3.set_xticks([-125, -124, -123, -122])
    ax3.set_xticklabels([-125, -124, -123, -122])
    ax3.text(.05,.05,'Volume Flux $[m^{3}s^{-1}]$', size=.8*fs, transform=ax3.transAxes,
        weight='bold', color=cq)

    pfun.add_coast(ax4, color='gray')
    ax4.axis(aaP)
    pfun.dar(ax4)
    ax4.set_yticklabels([])
    ax4.tick_params(labelsize=.8*fs)
    ax4.text(.95,.9,'(d)', size=fs, transform=ax4.transAxes, ha='right', weight='bold')
    ax4.set_xlabel('Longitude', size=.8*fs)
    ax4.set_xticks([-123, -122.5])
    ax4.set_xticklabels([-123, -122.5])
        
    fig.tight_layout()
    plt.savefig(out_fn)
    
plt.show()


    

    