"""
This code loads in the results of a tracker release for which we have ROMS
diagnostics and averages.  It collects the u and v momentum budget terms
along selected particle paths and saves them for later analysis.

Performance: 1 minute per day on my mac (20 particles)
           : 17 minutes per day on boiler (9000 particles)

For testing run as:

run Ldyn_gather -testing True

"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun

Ldir = Lfun.Lstart()

import numpy as np
import netCDF4 as nc
from scipy.spatial import cKDTree
import pickle
import pandas as pd
from time import time
from datetime import datetime, timedelta

import Ldyn_functions as Ldf
from importlib import reload
reload(Ldf)

# command line inputs
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-exp_name', type=str, default='EJdF3d_3d_up4')
parser.add_argument('-release', type=str, default='2018.05.15')
parser.add_argument('-testing', type=zfun.boolean_string, default=False)
parser.add_argument('-verbose', type=zfun.boolean_string, default=False)
parser.add_argument('-roms', type=str, default='roms3')
args = parser.parse_args()
exp_name = args.exp_name
testing = args.testing
verbose = args.verbose
rds = args.release

t_dir = Ldir['LOo'] + 'tracks/' + exp_name + '/'
EI = Lfun.csv_to_dict(t_dir + 'exp_info.csv')

Ldir['roms'] = Ldir[args.roms]

t_fn = 'release_' + rds + '.nc'

out_dir = t_dir + 'Ldyn_' + rds + '/'
if verbose:
    print(t_fn)
    print(rds)
    print(out_dir)
    
Lfun.make_dir(out_dir, clean = True)

if testing == False:
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
    t_ds.close()

    # Choose the "winners"
    if 'EJdF3d' in exp_name:
        seg_list = Ldf.seg_dict['PSTrim']
    
    plon = lon[-1,:]
    plat = lat[-1,:]
    G = zrfun.get_basic_info(EI['fn00'], only_G=True)
    glon = zfun.fillit(G['lon_rho'])
    glat = zfun.fillit(G['lat_rho'])
    mask, imask = Ldf.get_imask(Ldir, seg_list, plon, plat, glon, glat)
    NPM = len(imask)

    # load pre-made trees
    tree_dir = Ldir['LOo'] + 'tracker_trees/' + EI['gridname'] + '/'
    # 2D
    xyT_rho = pickle.load(open(tree_dir + 'xyT_rho.p', 'rb'))
    xyT_u = pickle.load(open(tree_dir + 'xyT_u.p', 'rb'))
    xyT_v = pickle.load(open(tree_dir + 'xyT_v.p', 'rb'))
    xyT_rho_un = pickle.load(open(tree_dir + 'xyT_rho_un.p', 'rb'))
    # 3D
    xyzT_rho = pickle.load(open(tree_dir + 'xyzT_rho.p', 'rb'))
    xyzT_u = pickle.load(open(tree_dir + 'xyzT_u.p', 'rb'))
    xyzT_v = pickle.load(open(tree_dir + 'xyzT_v.p', 'rb'))
    xyzT_w = pickle.load(open(tree_dir + 'xyzT_w.p', 'rb'))

    # ROMS file stuff

    # make lists of all the associated his, dia and avg files for this time span
    Ldir['gtagex'] = EI['gtagex']
    date_string0 = EI['ds_first_day']
    dt0 = datetime.strptime(date_string0, '%Y.%m.%d')
    dt1 = dt0 + timedelta(days=int(EI['days_to_track'])-1)
    date_string1 = dt1.strftime('%Y.%m.%d')
    his_list = Lfun.get_fn_list('hourly', Ldir, date_string0, date_string1)

    date_list = Lfun.date_list_utility(dt0, dt1)
    roms_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
    avg_list = []
    dia_list = []
    date_list = Lfun.date_list_utility(dt0, dt1)
    for dl in date_list:
        f_string = 'f' + dl
        for nhis in range(1,25):
            nhiss = ('0000' + str(nhis))[-4:]
            avg_list.append(roms_dir + f_string + '/ocean_avg_' + nhiss + '.nc')
            dia_list.append(roms_dir + f_string + '/ocean_dia_' + nhiss + '.nc')

    # also make lists of datetimes
    his_dt_list = []
    for fn in his_list:
        T = zrfun.get_basic_info(fn, only_T=True)
        his_dt_list.append(T['tm'])
    avg_dt_list = []
    for fn in avg_list:
        T = zrfun.get_basic_info(fn, only_T=True)
        avg_dt_list.append(T['tm'])
    dia_dt_list = []
    for fn in dia_list:
        T = zrfun.get_basic_info(fn, only_T=True)
        dia_dt_list.append(T['tm'])

    # things for nearest neighbor interpolation
    Ldir = Lfun.Lstart()
    G, S, T = zrfun.get_basic_info(EI['fn00'])
    #Maskr = G['mask_rho'] # True over water
    Masku = G['mask_u'] # True over water
    Maskv = G['mask_v'] # True over water
    #Maskr3 = np.tile(G['mask_rho'].reshape(1,G['M'],G['L']),[S['N'],1,1])
    Masku3 = np.tile(G['mask_u'].reshape(1,G['M'],G['L']-1),[S['N'],1,1])
    Maskv3 = np.tile(G['mask_v'].reshape(1,G['M']-1,G['L']),[S['N'],1,1])
    #Maskw3 = np.tile(G['mask_rho'].reshape(1,G['M'],G['L']),[S['N']+1,1,1])

    # now go through each track and get diagnostics

    # first average the winners onto the half hour
    wlon = (lon[:-1,imask] + lon[1:,imask])/2
    wlat = (lat[:-1,imask] + lat[1:,imask])/2
    wcs = (cs[:-1,imask] + cs[1:,imask])/2

    if testing == True:
        dia_list = dia_list[:2]
        dia_dt_list = dia_dt_list[:2]

    a = pd.DataFrame(index=dia_dt_list, columns=imask)
    u_dict = {}
    v_dict = {}
    a_list = ['accel','hadv', 'vadv','cor','prsgrd','prsgrd0','vvisc']
    u_list = ['u_'+item for item in a_list]
    v_list = ['v_'+item for item in a_list]
    for vn in u_list:
        u_dict[vn] = a.copy()
    for vn in v_list:
        v_dict[vn] = a.copy()

    tt0 = time()
    ii = 0
    for fn in dia_list:
        print('Working on %s' % (dia_dt_list[ii].strftime('%Y.%m.%d %H:%M')))
        sys.stdout.flush()
        ds = nc.Dataset(fn)
        for vn in u_list:
            if 'prsgrd0' in vn:
                A = ds[vn.replace('0','')][0,:,:,:]
                Af = A[Masku3].data
                # also get the pressure gradient at the surface
                xys = np.array((wlon[ii,:],wlat[ii,:],0*wcs[ii,:])).T
            else:
                A = ds[vn][0,:,:,:]
                Af = A[Masku3].data
                xys = np.array((wlon[ii,:],wlat[ii,:],wcs[ii,:])).T
            Ai = Af[xyzT_u.query(xys, n_jobs=-1)[1]]
            u_dict[vn].loc[dia_dt_list[ii],:] = Ai
            
        for vn in v_list:
            if 'prsgrd0' in vn:
                A = ds[vn.replace('0','')][0,:,:,:]
                Af = A[Maskv3].data
                # also get the pressure gradient at the surface
                xys = np.array((wlon[ii,:],wlat[ii,:],0*wcs[ii,:])).T
            else:
                A = ds[vn][0,:,:,:]
                Af = A[Maskv3].data
                xys = np.array((wlon[ii,:],wlat[ii,:],wcs[ii,:])).T
            Ai = Af[xyzT_v.query(xys, n_jobs=-1)[1]]
            v_dict[vn].loc[dia_dt_list[ii],:] = Ai
        ii += 1
    print('Filling u_dict and v_dict took %0.1f sec' % (time()-tt0))


    # save results
    pickle.dump(u_dict, open(out_dir + 'Ldyn_u_dict.p', 'wb'))
    pickle.dump(v_dict, open(out_dir + 'Ldyn_v_dict.p', 'wb'))
    pickle.dump(imask, open(out_dir + 'Ldyn_imask.p', 'wb'))

    if False:
        # testing balance: RESULT it looks perfect
        print('\n** x-mom terms **')
        x = 0
        for vn in u_list:
            aa = u_dict[vn].iloc[0,0]
            print('  %s = %0.2e' % (vn, aa))
            if 'accel' in vn:
                x += aa
            elif 'prsgrd0' in vn:
                pass
            else:
                x -= aa
        print(' sum = %0.2e' % (x))
    
        print('\n** y-mom terms **')
        y = 0
        for vn in v_list:
            aa = v_dict[vn].iloc[0,0]
            print('  %s = %0.2e' % (vn, aa))
            if 'accel' in vn:
                y += aa
            elif 'prsgrd0' in vn:
                pass
            else:
                y -= aa
        print(' sum = %0.2e' % (x))

    
