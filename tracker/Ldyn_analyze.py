"""
This code loads in the results of Ldyn_gather and analyzes the
along-track momentum balance of the "winners."

Performance:

For testing run as:

"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun

sys.path.append(os.path.abspath('../plotting'))
import pfun

Ldir = Lfun.Lstart()

import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
from time import time
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import Ldyn_functions as Ldf
from importlib import reload
reload(Ldf)

# command line inputs
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-exp_name', type=str, default='EJdF3d_3d_up4')
parser.add_argument('-testing', type=zfun.boolean_string, default=False)
args = parser.parse_args()
exp_name = args.exp_name
testing = args.testing

t_dir = Ldir['LOo'] + 'tracks/' + exp_name + '/'
EI = Lfun.csv_to_dict(t_dir + 'exp_info.csv')
p_dir = t_dir + 'Lplots/'
Lfun.make_dir(p_dir, clean=True)

rds_list = ['2018.05.15', '2018.05.22', '2018.05.29', '2018.06.05']
if testing:
    rds_list = [rds_list[0]]
    
for rds in rds_list:
    
    # out name for plot
    pname = p_dir + 'Ldyn_' + rds + '.png'

    t_fn = 'release_' + rds + '.nc'
    # get the tracker output for this release
    t_ds = nc.Dataset(t_dir + t_fn)
    NT, NP = t_ds['lon'].shape
    # get a list of datetimes
    ot_vec = t_ds['ot'][:]
    t_dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]
    # Gather particle data; packed [time, particle #]
    lon = t_ds['lon'][:]
    lat = t_ds['lat'][:]
    cs = t_ds['cs'][:]
    u = t_ds['u'][:]
    v = t_ds['v'][:]
    t_ds.close()

    # load results of Ldyn_gather.py
    
    u_dict = pickle.load(open(t_dir + 'Ldyn_' + rds +'/Ldyn_u_dict.p', 'rb'))
    v_dict = pickle.load(open(t_dir + 'Ldyn_' + rds +'/Ldyn_v_dict.p', 'rb'))
    imask = pickle.load(open(t_dir + 'Ldyn_' + rds +'/Ldyn_imask.p', 'rb'))

    # base names of the various acceleration diagnostics
    # (prsgrd0 is prsgrd at the surface)
    a_list = ['accel','hadv', 'vadv','cor','prsgrd','prsgrd0','vvisc']

    if testing == True:
        imask = imask[::10]
        for vn in a_list:
            u_dict['u_'+vn] = u_dict['u_'+vn][imask]
            v_dict['v_'+vn] = v_dict['v_'+vn][imask]

    # make sure everything is numeric
    for vn in a_list:
        for cn in u_dict['u_'+vn].columns:
            u_dict['u_'+vn][cn] = pd.to_numeric(u_dict['u_'+vn][cn])
        for cn in v_dict['v_'+vn].columns:
            v_dict['v_'+vn][cn] = pd.to_numeric(v_dict['v_'+vn][cn])

    NTW, NPW = u_dict['u_accel'].shape

    # find x, y, etc. for the particles
    wlon = zfun.fillit(lon[:,imask])
    wlat = zfun.fillit(lat[:,imask])
    wcs = zfun.fillit(cs[:,imask])
    wu = zfun.fillit(u[:,imask])
    wv = zfun.fillit(v[:,imask])
    x, y = zfun.ll2xy(wlon,wlat, wlon.mean(), wlat.mean())

    # rotate into along-mean-track momeuntum budget
    # meaning there is one direction for each track,
    # defined by its end points
    theta = np.zeros((NTW, NPW))
    for ii in range(x.shape[1]):
        DX = x[-1,ii] - x[0,ii]
        DY = y[-1,ii] - y[0,ii]
        theta[:,ii] = np.arctan2(DY,DX)
    sin_th = np.sin(theta)
    cos_th = np.cos(theta)

    # interpolate velocities onto the half hour (where we have diagnostics)
    # and rotate to track directions
    uu = (wu[1:,:]+wu[:-1,:])/2
    vv = (wv[1:,:]+wv[:-1,:])/2
    u_r = cos_th*uu + sin_th*vv
    v_r = cos_th*vv - sin_th*uu

    udf = pd.DataFrame(u_r, index=u_dict['u_accel'].index, columns=u_dict['u_accel'].columns)
    vdf = pd.DataFrame(v_r, index=u_dict['u_accel'].index, columns=u_dict['u_accel'].columns)

    s_dict = {}
    n_dict = {}
    # "s" stand's for along-track, "n" stands for cross-track
    for vn in a_list:
        s_dict['s_'+vn] = cos_th*u_dict['u_'+vn] + sin_th*v_dict['v_'+vn]
        n_dict['n_'+vn] = cos_th*v_dict['v_'+vn] - sin_th*u_dict['u_'+vn]
    
    # form along-track running-temporal means
    si_dict = {}
    ni_dict = {}
    for vn in a_list:
        si_dict['s_'+vn] = (s_dict['s_'+vn])#.cumsum(axis=0) / np.ones((NTW, NPW)).cumsum(axis=0)
        ni_dict['n_'+vn] = (n_dict['n_'+vn])#.cumsum(axis=0) / np.ones((NTW, NPW)).cumsum(axis=0)
    
    # form more informative groupings of terms
    SI_dict = {}
    SI_dict['accel'] = si_dict['s_accel'] - si_dict['s_hadv'] - si_dict['s_vadv'] # LHS
    SI_dict['CAP'] = si_dict['s_prsgrd'] + si_dict['s_cor'] - (si_dict['s_accel'] - si_dict['s_hadv'] - si_dict['s_vadv']) # RHS
    SI_dict['CAP0'] = si_dict['s_prsgrd0'] + si_dict['s_cor'] - (si_dict['s_accel'] - si_dict['s_hadv'] - si_dict['s_vadv']) # RHS
    SI_dict['fric'] = si_dict['s_vvisc'] # RHS

    # filter
    for vn in SI_dict.keys():
        # also convert to m/s per day
        SI_dict[vn].loc[:,:] = zfun.filt_godin_mat(SI_dict[vn].to_numpy())*86400
    udf.loc[:,:] = zfun.filt_godin_mat(udf.to_numpy())
    vdf.loc[:,:] = zfun.filt_godin_mat(vdf.to_numpy())

    # PLOTTING
    plt.close('all')
    fs = 14
    plt.rc('font', size=fs)

    if False:
        # plot individual tracks, one per figure
        for ii in range(10):
    
            fig = plt.figure(figsize=(16,8))
    
            ax = fig.add_subplot(121)
            ax.plot(wlon[:,ii], wlat[:,ii],'-b', lw=3)
            ax.plot(wlon[0,ii], wlat[0,ii],'ok', ms=20)
            ax.plot(wlon[-1,ii], wlat[-1,ii],'*y', ms=20)
            pfun.add_coast(ax)
            pad = .1
            ax.axis([wlon[:,ii].min()-pad, wlon[:,ii].max()+pad, wlat[:,ii].min()-pad, wlat[:,ii].max()+pad])
            pfun.dar(ax)
    
            ax = fig.add_subplot(222)
            for vn in SI_dict.keys():
                SI_dict[vn][imask[ii]].plot(ax=ax, legend=True, label=vn, grid=True)
            ax.set_xticklabels([])
            ax.set_xticklabels([], minor=True)
        
            ax = fig.add_subplot(224)
            udf[imask[ii]].plot(ax=ax, legend=True, label='Along-track Velocity', grid=True)
    else:
        # plot many tracks on a figure
        fig = plt.figure(figsize=(16,8))
    
        al = .1
    
        # map
        ax = fig.add_subplot(121)
        ax.plot(wlon[:,:], wlat[:,:],'-g', lw=1, alpha=al)
        ax.plot(wlon[0,:], wlat[0,:],'o',c='orange', ms=4, alpha=.8)
        ax.plot(wlon[-1,:], wlat[-1,:],'o', c='m', ms=5, alpha=.8)
        pfun.add_coast(ax)
        pad = .1
        #ax.axis([wlon.min()-pad, wlon.max()+pad, wlat.min()-pad, wlat.max()+pad])
        ax.axis([-124, -122, 47, 49])
        ax.set_xticks([-124, -123, -122])
        ax.set_yticks([47, 48, 49])
        pfun.dar(ax)
        ax.text(.05,.2,'START', c='orange', transform=ax.transAxes, weight='bold')
        ax.text(.05,.1,'FINISH', c='m', transform=ax.transAxes, weight='bold')
        ax.set_title(exp_name + ' ' + rds)

        # forces/unit mass
        ax = fig.add_subplot(222)
        SI_dict['fric'].plot(ax=ax, legend=False, grid=True, alpha=al, c='k', lw=.5)
        SI_dict['fric'].mean(axis=1).plot(ax=ax, legend=False, grid=True, alpha=1, lw=3, c='r')
        SI_dict['CAP'].mean(axis=1).plot(ax=ax, legend=False, grid=True, alpha=1, lw=3, c='b')
        SI_dict['CAP0'].mean(axis=1).plot(ax=ax, legend=False, grid=True, alpha=1, lw=3, c='c')
        ax.set_xticklabels([])
        ax.set_xticklabels([], minor=True)
        ax.text(.05,.1,r'Friction $[m\ s^{-1}\ day^{-1}]$',transform=ax.transAxes, c='r')
        ax.text(.05,.2,r'$CAP\ [m\ s^{-1}\ day^{-1}]$',transform=ax.transAxes, c='b')
        ax.text(.05,.3,r'$CAP_{0}\ [m\ s^{-1}\ day^{-1}]$',transform=ax.transAxes, c='c')
        ax.axhline(c='k', ls='--')
        ax.set_ylim(-1, 1)
    
        # velocity
        ax = fig.add_subplot(224)
        udf.plot(ax=ax, legend=False, grid=True, alpha=al, c='k', lw=.5)
        udf.mean(axis=1).plot(ax=ax, legend=False, grid=True, alpha=1, lw=3, c='k')
        ax.text(.05,.1,r'Along-track Velocity $[m\ s^{-1}]$',transform=ax.transAxes)
        ax.axhline(c='k', ls='--')
        ax.set_ylim(-.3, .3)
    
    
    plt.show()
    fig.savefig(pname)
    plt.rcdefaults()

