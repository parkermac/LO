"""
Plot results of a particle tracking experiment, specific to experiments about
why water is exchanged across the Admiralty Inlet sill.
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun
import zfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import seawater as sw

Ldir = Lfun.Lstart()

# command line inputs
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-exp_name', type=str, default='EJdF3d_3d_up4')
#parser.add_argument('-release', type=str, default='2018.05.15')
args = parser.parse_args()
exp_name = args.exp_name

t_dir = Ldir['LOo'] + 'tracks/' + exp_name + '/'
EI = Lfun.csv_to_dict(t_dir + 'exp_info.csv')
p_dir = t_dir + 'plots/'
Lfun.make_dir(p_dir, clean=True)

for rds in ['2018.05.15', '2018.05.22', '2018.05.29', '2018.06.05']:

    t_fn = 'release_' + rds + '.nc'
    
    # place for plots
    out_name = p_dir + t_fn.replace('.nc','.png')
    
    dsr = nc.Dataset(t_dir + t_fn)

    NT, NP = dsr['lon'].shape

    # get a list of datetimes
    ot_vec = dsr['ot'][:]
    dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]
    t = (ot_vec - ot_vec[0])/3600

    # Gather particle data
    # packed [time, particle #]
    lon = dsr['lon'][:]
    lat = dsr['lat'][:]
    z = dsr['z'][:]
    h = dsr['h'][:]
    u = dsr['u'][:]
    v = dsr['v'][:]
    w = dsr['w'][:]
    salt = dsr['salt'][:]
    temp = dsr['temp'][:]
    cs = dsr['cs'][:]
    zeta = dsr['zeta'][:]

    z = cs*h

    dsr.close()

    # Choose the "winners"
    import Ldyn_functions as Ldf
    from importlib import reload
    reload(Ldf)
    # Choose the "winners"
    if 'EJdF3d' in exp_name:
        seg_list = Ldf.seg_dict['PSTrim']
    plon = lon[-1,:]
    plat = lat[-1,:]
    dsg = nc.Dataset(t_dir + 'grid.nc')
    glon = zfun.fillit(dsg['lon_rho'][:])
    glat = zfun.fillit(dsg['lat_rho'][:])
    mask, imask = Ldf.get_imask(Ldir, seg_list, plon, plat, glon, glat)
    NPM = len(imask)

    # PLOTTING
    #plt.close('all')
    fs = 16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(18,10))

    # Map
    ax = fig.add_subplot(121)
    ax.plot(lon[0,:], lat[0,:], '.k', alpha=.1)
    ax.plot(lon[0,imask], lat[0,imask], '.c')
    ax.plot(lon[-1,imask], lat[-1,imask], '.b')
    pfun.dar(ax)
    pfun.add_coast(ax)
    aa = [-124, -122, 47, 49]
    ax.axis(aa)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-124, -123, -122])
    ax.set_yticks([47, 48, 49])
    ax.set_title(exp_name + ' ' + rds)
    ax.text(.05, .1, 'Winners = %0d%%' % (int(100*NPM/NP)), weight='bold',
        transform=ax.transAxes)

    # Histogram
    ax = fig.add_subplot(222)
    #zall = z.flatten()
    z0 = z[0,mask]
    z1 = z[-1,mask]
    bins=np.linspace(-300,0,30 + 1)
    #
    counts, obins = np.histogram(z0, bins=bins)
    ax.hist(bins[:-1], bins, weights=counts/NPM,
        orientation='horizontal', rwidth=.7, color='c')
    #
    counts, obins = np.histogram(z1, bins=bins)
    ax.hist(bins[:-1], bins, weights=counts/NPM,
        orientation='horizontal', rwidth=.4, color='b')
    ax.set_xlabel('Fraction')
    ax.set_ylabel('Z [m]')

    # Theta-S plot
    ax = fig.add_subplot(224)
    ax.plot(salt[0,mask], temp[0,mask],'.c', alpha=.7)
    ax.plot(salt[-1,mask], temp[-1,mask],'+b', alpha=.4)
    ax.set_xlabel('Salinity [g/kg]')
    ax.set_ylabel('Potential Temp. [degC]')
    # overlay potential density contours
    aa = [24, 34, 4, 18]
    ss = np.linspace(aa[0],aa[1],100)
    tt = np.linspace(aa[2],aa[3],100)
    SS, TT = np.meshgrid(ss,tt)
    RR = sw.dens0(SS,TT)
    CS = ax.contour(SS, TT, RR-1000)
    ax.clabel(CS, fmt='%d')
    ax.axis(aa)

    #plt.show()
    fig.tight_layout()
    fig.savefig(out_name)
    
    plt.rcdefaults()


