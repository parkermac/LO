"""
Calculate incident and reflected wave amplitudes on all sections.

RESULT: this appears to give unrelaible results, although I don't know why.
The way I test if a results is "correct" is to compare the net tidal energy
flux due to the incident and reflected waves (FF_alt) to the original flux (F)
and the constituent-reconstructed flux (FF).  In general F and FF compare
well regardless of the sign I use for Eg and Ug.  The FF_alt is also good
(and identical to FF) but ONLY if I use a real value of a.  For any complex
a FF_alt != FF.

Curious: If I use real a, the net is always good, but the incident
and reflected parts vary with a (but not their sum).  Is a = np.sqrt(g/H)
more "correct?"

This also makes a nice map, if I ever get the answer to be satisfactory.

"""

from pathlib import Path
import sys
pth = Path(__file__).absolute().parent.parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# imports
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np
import cmath
import matplotlib.pyplot as plt
import pytide
from datetime import datetime

import Lfun
import plotting_functions as pfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun

# set input directory
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
year = 2018
dates = str(year) + '.01.01_' + str(year) + '.12.31'
in_dir = in_dir0 / ('harmonics_' + dates)

# get section info
sect_df = tef_fun.get_sect_df(gridname)

# loop over all sections
sect_list = list(sect_df.index)

testing = True
if testing:
    sect_list = ['ai1']

g = 9.8
rho = 1025
hm_e_dict = pickle.load(open(in_dir / 'hm_e_dict.p', 'rb'))
hm_u_dict = pickle.load(open(in_dir / 'hm_u_dict.p', 'rb'))
F_df = pd.DataFrame(index=sect_list)
for sect_name in sect_list:    
    hm_e = hm_e_dict[sect_name]
    hm_u = hm_u_dict[sect_name]
    H = hm_e['H']
    A0 = hm_e['A0']
    F = hm_e['F']
    
    # get section info
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    lon = (x0 + x1)/2
    lat = (y0 + y1)/2
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
    
    # net tidal energy flux
    FFp = 0
    FFm = 0
    FF = 0
    FF_alt = 0
    clist = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
    #clist = ['M2']

    # Use the pytide module to get the same information.
    # The tutorial is helpful:
    # https://pangeo-pytide.readthedocs.io/en/latest/tutorial.html
    wt = pytide.WaveTable(clist)
    f, vu = wt.compute_nodal_modulations([datetime(year,1,1)])
    
    for cons in clist:
        E = hm_e.A[hm_e.name == cons][0]
        U = hm_u.A[hm_u.name == cons][0]
        Eg = np.pi * hm_e.g[hm_e.name == cons][0] / 180
        Ug = np.pi * hm_u.g[hm_u.name == cons][0] / 180
        
        this_f = f[clist.index(cons)]
        this_vu = vu[clist.index(cons)]
        
        EE = cmath.rect(this_f*E, -Eg)
        UU = cmath.rect(this_f*U, -Ug)
        
        Fc = 0.5*rho*g*A0 * (EE.real*UU.real + EE.imag*UU.imag)
        
        a = np.sqrt(g/H)/np.sqrt(1 + 1j) # use 1j for R/om = 1
        
        Ep = (EE + UU/a)/2
        Em = (EE - UU/a)/2
        Up =  Ep*a
        Um =  Em*a
        Fp = 0.5*rho*g*A0 * (Ep.real*Up.real + Ep.imag*Up.imag)
        Fm = - 0.5*rho*g*A0 * (Em.real*Um.real + Em.imag*Um.imag)
        if testing:
            print('%s: Fp = %8.1f, Fm = %8.1f [MW]' % (cons, Fp/1e6, Fm/1e6))
        
        Fc_alt = Fp + Fm
        FFp += Fp
        FFm += Fm
        FF += Fc
        FF_alt += Fc_alt
        # FF should match FF_alt
            
    F_df.loc[sect_name, 'FFp'] = FFp/1e6
    F_df.loc[sect_name, 'FFm'] = FFm/1e6
    F_df.loc[sect_name, 'FF'] = FF/1e6
    F_df.loc[sect_name, 'F'] = F/1e6
    F_df.loc[sect_name, 'lon'] = lon
    F_df.loc[sect_name, 'lat'] = lat
    F_df.loc[sect_name, 'sdir'] = sdir
    
    if testing:
        print('\n%s: F = %0.1f, FF = %0.1f, FF_alt = %0.1f [MW]' %
            (sect_name, F/1e6, FF/1e6, FF_alt/1e6))

if not testing:
    
    # PLOTTING
    plt.close('all')
    fs = 16 # fontsize
    ffs = .7*fs # flux text fontsize
    alpha = .5 # transparency for arrows
    cf = 'purple' # color for incident tidal energy flux
    cq = 'g' # color for reflected tidal energy flux

    pfun.start_plot(fs=fs, figsize=(18,12))
    fig = plt.figure()

    # Angles by which to adjust the plotting of flux at each section.
    # They are defined -90:90 degrees, positive counter-clockwise, to define the
    # direction of an arrow pointing along the local thalweg, relative to the
    # normal to the section (positive East or North).
    th_dict = {
        'jdf1': -15, 'jdf2': -15, 'jdf3': -5, 'jdf4': 0,
        'sog1': 45, 'sog2': 55, 'sog3': -45, 'sog4': 15, 'sog5': 0,
        'sji1': 0, 'sji2': 20,
        'dp': 0,
        'ai1': -45, 'ai2': 30, 'ai3': 0, 'ai4': -30,
        'wb1': -30, 'wb2': 45, 'wb3': 0, 'wb4': 15,
        'hc1': 10, 'hc2': -45, 'hc3': 10, 'hc4': -50, 'hc5': -45, 'hc6': -30, 'hc7': 0, 'hc8': 25,
        'mb1': 0, 'mb2': -10, 'mb3': 0, 'mb4': 10, 'mb5': -45,
        'tn1': 0, 'tn2': 0, 'tn3': -15,
        'ss1': 45, 'ss2': -30, 'ss3': 0,
    }

    # # axis limits
    x0 = -125.5; x1 = -122; y0 = 47; y1 = 50.5 # Salish Sea
    x00 = -123.3; x11 = -122.2; y00 = 47; y11 = 48.5 # Puget Sound
    aaS = [x0, x1, y0, y1]
    aaP = [x00, x11, y00, y11]

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    #
    for sn in F_df.index:
        sx = F_df.loc[sn,'lon']
        sy = F_df.loc[sn,'lat']
        sdir = F_df.loc[sn, 'sdir']
    
        FF = F_df.loc[sn,'FF']
        FFp = F_df.loc[sn,'FFp']
        FFm = F_df.loc[sn,'FFm']
    
        if FF < 0:
            FFp, FFm = (FFm, FFp)
        
        if sdir == 'EW':
            om = 90 + th_dict[sn]
        elif sdir == 'NS':
            om = th_dict[sn]
        sth = np.sin(np.pi*om/180)
        cth = np.cos(np.pi*om/180)
    
        vFF0 = cth*FF
        vFF1 = sth*FF
    
        vFFp0 = cth*FFp
        vFFp1 = sth*FFp
    
        vFFm0 = cth*FFm
        vFFm1 = sth*FFm

        scl1 = 50000 # a vector of length 1 will be 1/scl of the length of the y-axis
        ax1.quiver(sx,sy, vFF0, vFF1, scale=scl1, scale_units='height', linewidths=1, color='k', alpha=1)
        ax1.quiver(sx,sy, vFFp0, vFFp1, scale=scl1, scale_units='height', linewidths=1, color=cf, alpha=alpha)
        ax1.quiver(sx,sy, vFFm0, vFFm1, scale=scl1, scale_units='height', linewidths=1, color=cq, alpha=alpha)
        ax1.plot(sx,sy,'ok')
        # if sx < x00 or sy > y11:
        #     ax1.text(sx, sy+.04, str(int(np.abs(F))), ha='center', va='center', size=ffs,
        #         weight='bold', color=cf, alpha=1)

        if sx > x00 and sy < y11:
            scl2 = 10000 # a vector of length 1 will be 1/scl of the length of the y-axis
            ax2.quiver(sx,sy, vFF0, vFF1, scale=scl2, scale_units='height', linewidths=1, color='k', alpha=1)
            ax2.quiver(sx,sy, vFFp0, vFFp1, scale=scl2, scale_units='height', linewidths=1, color=cf, alpha=alpha)
            ax2.quiver(sx,sy, vFFm0, vFFm1, scale=scl2, scale_units='height', linewidths=1, color=cq, alpha=alpha)
            ax2.plot(sx,sy,'ok')
            # ax2.text(sx, sy+.02, str(int(np.abs(F))), ha='center', va='center', size=ffs,
            #     weight='bold', color=cf, alpha=1)

    pfun.add_coast(ax1, color='gray')
    ax1.axis(aaS)
    pfun.dar(ax1)
    ax1.tick_params(labelsize=.8*fs)
    ax1.text(.95,.9,'(a)', size=fs, transform=ax1.transAxes, ha='right', weight='bold')
    ax1.set_xlabel('Longitude', size=.8*fs)
    ax1.set_ylabel('Latitude', size=.8*fs)
    ax1.set_xticks([-125, -124, -123, -122])
    ax1.text(.05,.12,'Tidal energy Flux $[MW]$', size=fs, transform=ax1.transAxes,
        weight='bold', color=cf)
    # legend
    sx = -124.5; sy = 47.15
    vFF0 = 2000; vFF1 = 0
    vFFp0 = 7000; vFFp1 = 0
    vFFm0 = -5000; vFFm1 = 0
    ax1.quiver(sx,sy, vFF0, vFF1, scale=scl1, scale_units='height', linewidths=1, color='k', alpha=1)
    ax1.quiver(sx,sy, vFFp0, vFFp1, scale=scl1, scale_units='height', linewidths=1, color=cf, alpha=alpha)
    ax1.quiver(sx,sy, vFFm0, vFFm1, scale=scl1, scale_units='height', linewidths=1, color=cq, alpha=alpha)
    ax1.plot(sx,sy,'ok')
    ax1.text(sx+.65, sy+.04, str(vFFp0), c=cf, weight='bold', alpha=alpha, ha='center')
    ax1.text(sx+.65, sy-.04, 'Incident', va='top', c=cf, weight='bold', alpha=alpha, ha='center')
    ax1.text(sx-.5, sy-.04, 'Reflected', va='top', c=cq, weight='bold', alpha=alpha, ha='center')
    ax1.text(sx+.1, sy-.04, 'Net', va='top', c='k', weight='bold', alpha=alpha, ha='center')

    
    pfun.add_coast(ax2, color='gray')
    ax2.axis(aaP)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude', size=.8*fs)
    ax2.set_yticks([47, 48])
    ax2.tick_params(labelsize=.8*fs)
    ax2.text(.95,.9,'(b)', size=fs, transform=ax2.transAxes, ha='right', weight='bold')
    ax2.set_xticks([-123, -122.5])
    # legend
    sx = -123; sy = 47.05
    vFF0 = 400; vFF1 = 0
    vFFp0 = 1200; vFFp1 = 0
    vFFm0 = -800; vFFm1 = 0
    ax2.quiver(sx,sy, vFF0, vFF1, scale=scl2, scale_units='height', linewidths=1, color='k', alpha=1)
    ax2.quiver(sx,sy, vFFp0, vFFp1, scale=scl2, scale_units='height', linewidths=1, color=cf, alpha=alpha)
    ax2.quiver(sx,sy, vFFm0, vFFm1, scale=scl2, scale_units='height', linewidths=1, color=cq, alpha=alpha)
    ax2.plot(sx,sy,'ok')
    ax2.text(sx+.2, sy+.02, str(vFFp0), c=cf, weight='bold', alpha=alpha)


    fig.suptitle(Ldir['gtagex'] + ': ' + in_dir.name.replace('harmonics_',''))
    fig.tight_layout()

    #fig.savefig(out_fn)
    plt.show()
    pfun.end_plot()
