"""
Plot bulk fluxes as a time series.

To test on mac:
run bulk_plot -gtx cas7_trapsV00_meV00 -ctag c0 -0 2017.01.01 -1 2017.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('bulk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])

if Ldir['testing'] == False:
    Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name.replace('.nc','') for item in in_dir.glob('*.nc')]

if Ldir['testing']:
    sect_list = ['ai7']
    from importlib import reload
    reload(tef_fun)
    
# grid info
g = xr.open_dataset(Ldir['grid'] / 'grid.nc')
h = g.h.values
h[g.mask_rho.values==0] = np.nan
xrho = g.lon_rho.values
yrho = g.lat_rho.values
xp, yp = pfun.get_plon_plat(xrho,yrho)
xu = g.lon_u.values
yu = g.lat_u.values
xv = g.lon_v.values
yv = g.lat_v.values

# PLOTTING
fs = 12
plt.close('all')
if Ldir['testing']:
    figsize = (12,8)
else:
    figsize = (21,12)
pfun.start_plot(fs=fs, figsize=figsize)

def add_qprism(ax):
    # add Qprism
    axqp = ax.twinx()
    axqp.plot(ot,tef_df['Q_prism'].to_numpy(),'-',
        color='c',linewidth=3,alpha=.4)
    axqp.text(.95,.9,r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$', color='c', 
        transform=ax.transAxes, ha='right',
        bbox=pfun.bbox)
    axqp.set_ylim(bottom=0)
    axqp.xaxis.label.set_color('c')
    axqp.tick_params(axis='y', colors='c')

for sect_name in sect_list:
    
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))

    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism'] = tef_df['qprism']/1000
                    
    # labels and colors
    ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
                'salt': r'Salinity $[g\ kg^{-1}]$',
                'QinDS': r'$Q_{in}\Delta S\ [10^{3}\ m^{3}s^{-1}\ g\ kg^{-1}]$'}
    p_color = 'r'
    m_color = 'b'
    lw = 2
        
    fig = plt.figure()
    
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=2) # Qin, Qout
    ax2 = plt.subplot2grid((3,3), (1,0), colspan=2) # Sin, Sout
    ax3 = plt.subplot2grid((3,3), (2,0), colspan=2) # Qin*DS and Qprism
    axmap = plt.subplot2grid((1,3), (0,2)) # map
    
    ot = bulk.time.values
    
    ax1.plot(ot,tef_df['Q_p'].to_numpy(), color=p_color, linewidth=lw)
    ax1.plot(ot,tef_df['Q_m'].to_numpy(), color=m_color, linewidth=lw)
    ax1.grid(True)    
    ax1.set_ylabel(ylab_dict['Q'])
    ax1.set_xlim(ot[0],ot[-1])
    
    qp = bulk['q'].values/1000
    qp[qp<0] = np.nan
    qm = bulk['q'].values/1000
    qm[qm>0]=np.nan
    sp = bulk['salt'].values
    sp[np.isnan(qp)] = np.nan
    sm = bulk['salt'].values
    sm[np.isnan(qm)]=np.nan
    
    alpha=.3
    ax1.plot(ot,qp,'or',alpha=alpha)
    ax1.plot(ot,qm,'ob',alpha=alpha)
    
    ax2.plot(ot,sp,'or',alpha=alpha)
    ax2.plot(ot,sm,'ob',alpha=alpha)
    
    ax2.plot(ot,tef_df['salt_p'].to_numpy(), color=p_color, linewidth=lw)
    ax2.plot(ot,tef_df['salt_m'].to_numpy(), color=m_color, linewidth=lw)
    ax2.grid(True)
    ax2.set_ylabel(ylab_dict['salt'])
    ax2.set_xlim(ot[0],ot[-1])
    
    # Add Qin*DS, which requires figuring out which direction is "in".
    # We use the direction that has the higher mean salinity.
    sp_mean = tef_df['salt_p'].mean() # Ignores NaN's
    sm_mean = tef_df['salt_m'].mean()
    if sp_mean > sm_mean:
        DS = tef_df['salt_p'].to_numpy()-tef_df['salt_m'].to_numpy()
        Qin = tef_df['Q_p'].to_numpy()
    elif sp_mean <= sm_mean:
        DS = tef_df['salt_m'].to_numpy()-tef_df['salt_p'].to_numpy()
        Qin = -tef_df['Q_m'].to_numpy()
    ax3.plot(ot,DS*Qin, color='k', linewidth=lw)
    ax3.grid(True)
    ax3.set_ylabel(ylab_dict['QinDS'])
    ax3.set_xlim(ot[0],ot[-1])
    add_qprism(ax3)
    if sp_mean > sm_mean:
        ax3.text(.05,.05,r'$Q_{in}$ is positive', transform=ax3.transAxes)
    elif sp_mean <= sm_mean:
        ax3.text(.05,.05,r'$Q_{in}$ is negative', transform=ax3.transAxes)
    
    # map
    sn = sect_name.replace('.p','')
    sinfo = sect_df.loc[sect_df.sn==sn,:]
    i0 = sinfo.iloc[0,:].i
    j0 = sinfo.iloc[0,:].j
    uv0 = sinfo.iloc[0,:].uv
    i1 = sinfo.iloc[-1,:].i
    j1 = sinfo.iloc[-1,:].j
    uv1 = sinfo.iloc[-1,:].uv
    if uv0=='u':
        x0 = xu[j0,i0]
        y0 = yu[j0,i0]
    elif uv0=='v':
        x0 = xv[j0,i0]
        y0 = yv[j0,i0]
    if uv1=='u':
        x1 = xu[j1,i1]
        y1 = yu[j1,i1]
    elif uv1=='v':
        x1 = xv[j1,i1]
        y1 = yv[j1,i1]
    axmap.plot([x0,x1],[y0,y1],'-c', lw=3)
    axmap.plot(x0,y0,'og', ms=10)
    pfun.add_coast(axmap)
    pfun.dar(axmap)
    axmap.pcolormesh(xp, yp, -h, vmin=-100, vmax=100,
        cmap='jet', alpha=.4)
    
    dx = x1-x0; dy = y1-y0
    xmid = x0 + dx/2; ymid = y0 + dy/2
    pad = np.max((np.sqrt(dx**2 + dy**2)*2,.1))
    axmap.axis([x0-pad, x1+pad, y0-pad, y1+pad])
    axmap.set_xlabel('Longitude [deg]')
    axmap.set_ylabel('Latitude [deg]')
    axmap.set_title(sn)
    axmap.text(xmid - np.max((dy,.05))/6, ymid + np.max((dx,.05))/6, '+', fontweight='bold', c='r', fontsize=20,
        ha='center',va='center')
                
    # fig.tight_layout()
    
    if Ldir['testing']:
        plt.show()
    else:
        plt.savefig(out_dir / (sect_name.replace('.p','') + '.png'))
        plt.close()

pfun.end_plot()
