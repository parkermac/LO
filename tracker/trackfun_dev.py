"""
Functions for particle tracking.

Relies on the existence of experiment-specific info in
LO_output/tracks/exp_info.csv.

"""
# setup (assume path to alpha set by calling code)
from lo_tools import Lfun, zfun, zrfun

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
import pickle
from time import time
import sys

verbose = False

# this is the full list of tracers we want to find on the track and write to output
tracer_list_full = ['salt', 'temp']#, 'oxygen']
# we trim the list below so that it only includes tracers that are present
# in the history files

# criterion for deciding if particles are on land
maskr_crit = 0.5 # (maskr = 1 in water, 0 on land) [0.5 seems good]

# NEW CODE for nearest neighbor interpolation
Ldir = Lfun.Lstart()

TR0 = Lfun.csv_to_dict(Ldir['LOo'] / 'tracks' / 'exp_info.csv')
# Read this in as TR0 instead of TR so as not to confuse it with the TR
# that we pass to get_tracks() below.  The only difference is that in
# TR0 everything has been turned into a string, whereas in get_tracks()
# we assume everything in TR has its original dtype as set in tracker.py.

ds = xr.open_dataset(TR0['fn00'])
tracer_list = [vn for vn in tracer_list_full if vn in ds.data_vars]
ds.close()

print(tracer_list)

G, S, T = zrfun.get_basic_info(TR0['fn00'])
Maskr = G['mask_rho']==1 # True over water
Masku = G['mask_u']==1 # True over water
Maskv = G['mask_v']==1 # True over water
Maskr3 = np.tile(Maskr.reshape(1,G['M'],G['L']),[S['N'],1,1])
Masku3 = np.tile(Masku.reshape(1,G['M'],G['L']-1),[S['N'],1,1])
Maskv3 = np.tile(Maskv.reshape(1,G['M']-1,G['L']),[S['N'],1,1])
Maskw3 = np.tile(Maskr.reshape(1,G['M'],G['L']),[S['N']+1,1,1])
# load pre-made trees
tree_dir = Ldir['LOo'] / 'tracker_trees' / TR0['gridname']
# 2D
xyT_rho = pickle.load(open(tree_dir / 'xyT_rho.p', 'rb'))
xyT_u = pickle.load(open(tree_dir / 'xyT_u.p', 'rb'))
xyT_v = pickle.load(open(tree_dir / 'xyT_v.p', 'rb'))
xyT_rho_un = pickle.load(open(tree_dir / 'xyT_rho_un.p', 'rb'))
# 3D
xyzT_rho = pickle.load(open(tree_dir / 'xyzT_rho.p', 'rb'))
xyzT_u = pickle.load(open(tree_dir / 'xyzT_u.p', 'rb'))
xyzT_v = pickle.load(open(tree_dir / 'xyzT_v.p', 'rb'))
xyzT_w = pickle.load(open(tree_dir / 'xyzT_w.p', 'rb'))
# the "f" below refers to flattened, which is the result of passing
# a Boolean array like Maskr to an array.
lonrf = G['lon_rho'][Maskr]
latrf = G['lat_rho'][Maskr]

maskr = Maskr.flatten()
# minimum grid sizes in degrees (assumes plaid lon,lat grid?)
dxg = np.diff(G['lon_rho'][0,:]).min()
dyg = np.diff(G['lat_rho'][:,0]).min()

def get_tracks(fn_list, plon0, plat0, pcs0, TR, trim_loc=False):
    """
    This is the main function doing the particle tracking.
    """
    # unpack items needed from TR
    surface = not TR['3d']
    turb = TR['turb']
    ndiv = TR['ndiv']
    windage = TR['windage']
    
    # get time vector of history files
    NT = len(fn_list)
    NTS = TR['sph']*(NT-1) + 1
    rot = np.nan * np.ones(NT)
    counter = 0
    for fn in fn_list:
        ds = xr.open_dataset(fn, decode_times=False)
        rot[counter] = ds.ocean_time.values[0]
        counter += 1
        ds.close
    sys.stdout.flush()
    delta_t_his = rot[1] - rot[0] # seconds between saves
    delt = delta_t_his/ndiv # time step in seconds
    delt_save = delta_t_his/TR['sph']
    rot_save = np.linspace(rot[0], rot[-1], NTS)

    # these lists are used internally to get other variables as needed
    # (the variables must all be present in the history files)
    vn_list_vel = ['u','v','w']
    vn_list_zh = ['zeta','h']
    vn_list_wind = ['Uwind','Vwind']
    vn_list_other = tracer_list + vn_list_zh + vn_list_vel
    if windage > 0:
        vn_list_other = vn_list_other + vn_list_wind
    # plist_main is what ends up written to output
    plist_main = ['lon', 'lat', 'cs', 'ot', 'z'] + vn_list_other
    
    # debugging
    vconst = True # set v to a constant, u and w to zero
    jx = True # use Jilian's modification of the pcs step

    # Step through times.
    #
    for counter_his in range(len(fn_list)-1):
        
        if verbose:
            print('  >> hour = ' + str(counter_his))
            sys.stdout.flush()
                
        it0 = TR['sph']*counter_his
        
        # get Datasets
        fn0 = fn_list[counter_his]
        fn1 = fn_list[counter_his+1]
        
        ds0 = xr.open_dataset(fn0)
        ds1 = xr.open_dataset(fn1)
        
        # prepare the fields for nearest neighbor interpolation
        tt0 = time()
        if counter_his == 0:
            h = G['h']
            hf = h[Maskr]
            if surface == True:
                u0 = ds0['u'][0,-1,:,:].values
                u1 = ds1['u'][0,-1,:,:].values
                uf0 = u0[Masku]
                uf1 = u1[Masku]
                v0 = ds0['v'][0,-1,:,:].values
                v1 = ds1['v'][0,-1,:,:].values
                vf0 = v0[Maskv]
                vf1 = v1[Maskv]
                wf0 = 0; wf1 = 0
                
                trf0_dict = dict()
                trf1_dict = dict()
                for vn in tracer_list:
                    tr0 = ds0[vn][0,-1,:,:].values
                    tr1 = ds1[vn][0,-1,:,:].values
                    trf0_dict[vn] = tr0[Maskr]
                    trf1_dict[vn] = tr1[Maskr]
                
                if windage > 0:
                    Uwind0 = ds0['Uwind'][0,:,:].values
                    Uwind1 = ds1['Uwind'][0,:,:].values
                    Uwindf0 = Uwind0[Maskr]
                    Uwindf1 = Uwind1[Maskr]
                    Vwind0 = ds0['Vwind'][0,:,:].values
                    Vwind1 = ds1['Vwind'][0,:,:].values
                    Vwindf0 = Vwind0[Maskr]
                    Vwindf1 = Vwind1[Maskr]
            else:
                u0 = ds0['u'][0,:,:,:].values
                u1 = ds1['u'][0,:,:,:].values

                if vconst:
                    u0 = np.zeros(u0.shape)
                    u1 = u0.copy()

                uf0 = u0[Masku3]
                uf1 = u1[Masku3]
                v0 = ds0['v'][0,:,:,:].values
                v1 = ds1['v'][0,:,:,:].values

                if vconst:
                    v0 = 0.1 * np.ones(v0.shape)
                    v1 = v0.copy()

                vf0 = v0[Maskv3]
                vf1 = v1[Maskv3]
                w0 = ds0['w'][0,:,:,:].values
                w1 = ds1['w'][0,:,:,:].values

                if vconst:
                    w0 = np.zeros(w0.shape)
                    w1 = w0.copy()

                wf0 = w0[Maskw3]
                wf1 = w1[Maskw3]

                if vconst:
                    uf00 = uf0.copy()
                    vf00 = vf0.copy()
                    wf00 = wf0.copy()
                
                trf0_dict = dict()
                trf1_dict = dict()
                for vn in tracer_list:
                    tr0 = ds0[vn][0,:,:,:].values
                    tr1 = ds1[vn][0,:,:,:].values
                    trf0_dict[vn] = tr0[Maskr3]
                    trf1_dict[vn] = tr1[Maskr3]
                    
                if turb == True:
                    AKs0_temp = ds0['AKs'][0,:,:,:].values
                    # modify top and bottom AKs to be non-negligible
                    AKs0_temp[0,:,:] = AKs0_temp[1,:,:]
                    AKs0_temp[-1,:,:] = AKs0_temp[-2,:,:]
                    AKs0 = AKs0_temp.copy()
                    AKs0[1:-1,:,:] = 0.25*AKs0_temp[:-2,:,:] + 0.5*AKs0_temp[1:-1,:,:] + 0.25*AKs0_temp[2:,:,:]
                    AKsf0 = AKs0[Maskw3]
                    #
                    AKs1_temp = ds1['AKs'][0,:,:,:].values
                    # modify top and bottom AKs to be non-negligible
                    AKs1_temp[0,:,:] = AKs1_temp[1,:,:]
                    AKs1_temp[-1,:,:] = AKs1_temp[-2,:,:]
                    AKs1 = AKs1_temp.copy()
                    AKs1[1:-1,:,:] = 0.25*AKs1_temp[:-2,:,:] + 0.5*AKs1_temp[1:-1,:,:] + 0.25*AKs1_temp[2:,:,:]
                    AKsf1 = AKs1[Maskw3]
                    #
                    # New 2022.11.14 use time-varying dz
                    zeta0 = ds0['zeta'].values #jx
                    zeta1 = ds1['zeta'].values #jx
                    zw0 = zrfun.get_z(G['h'], zeta0, S, only_w=True) #jx
                    zw1 = zrfun.get_z(G['h'], zeta1, S, only_w=True) #jx
                    dz0 = np.diff(zw0, axis=0) #jx
                    dz1 = np.diff(zw1, axis=0) #jx
                    dKdz0 = np.diff(AKs0, axis=0)/dz0 #jx
                    dKdz1 = np.diff(AKs1, axis=0)/dz1 #jx
                    dKdzf0 = dKdz0[Maskr3]
                    dKdzf1 = dKdz1[Maskr3]
                    
                    
            z0 = ds0['zeta'][0,:,:].values
            z1 = ds1['zeta'][0,:,:].values
            zf0 = z0[Maskr]
            zf1 = z1[Maskr]
            if verbose:
                print('   > Prepare fields for tree %0.4f sec' % (time()-tt0))
        else:
            # subsequent time steps
            if surface == True:
                u1 = ds1['u'][0,-1,:,:].values
                uf0 = uf1.copy()
                uf1 = u1[Masku]
                v1 = ds1['v'][0,-1,:,:].values
                vf0 = vf1.copy()
                vf1 = v1[Maskv]
                wf0 = 0; wf1 = 0
                
                for vn in tracer_list:
                    tr1 = ds1[vn][0,-1,:,:].values
                    trf0_dict[vn] = trf1_dict[vn].copy()
                    trf1_dict[vn] = tr1[Maskr]
                
                if windage > 0:
                    Uwind1 = ds1['Uwind'][0,:,:].values
                    Uwindf0 = Uwindf1
                    Uwindf1 = Uwind1[Maskr]
                    Vwind1 = ds1['Vwind'][0,:,:].values
                    Vwindf0 = Vwindf1
                    Vwindf1 = Vwind1[Maskr]
            else:
                if TR['no_advection'] == False:
                    u1 = ds1['u'][0,:,:,:].values
                    uf0 = uf1.copy()
                    uf1 = u1[Masku3]
                    v1 = ds1['v'][0,:,:,:].values
                    vf0 = vf1.copy()
                    vf1 = v1[Maskv3]
                    w1 = ds1['w'][0,:,:,:].values
                    wf0 = wf1.copy()
                    wf1 = w1[Maskw3]

                    if vconst:
                        uf0 = uf00.copy()
                        uf1 = uf00.copy()
                        vf0 = vf00.copy()
                        vf1 = vf00.copy()
                        wf0 = wf00.copy()
                        wf1 = wf00.copy()
                    
                    for vn in tracer_list:
                        tr1 = ds1[vn][0,:,:,:].values
                        trf0_dict[vn] = trf1_dict[vn].copy()
                        trf1_dict[vn] = tr1[Maskr3]
                    
                if turb == True:
                    AKsf0 = AKsf1.copy()
                    #
                    AKs1_temp = ds1['AKs'][0,:,:,:].values
                    # modify top and bottom AKs to be non-negligible
                    AKs1_temp[0,:,:] = AKs1_temp[1,:,:]
                    AKs1_temp[-1,:,:] = AKs1_temp[-2,:,:]
                    AKs1 = AKs1_temp.copy()
                    AKs1[1:-1,:,:] = 0.25*AKs1_temp[:-2,:,:] + 0.5*AKs1_temp[1:-1,:,:] + 0.25*AKs1_temp[2:,:,:]
                    AKsf1 = AKs1[Maskw3]
                    #
                    dKdzf0 = dKdzf1.copy()
                    # new 2022.11.14 Use time-varying dz
                    zeta1 = ds1['zeta'].values #jx
                    zw1 = zrfun.get_z(G['h'], zeta1, S, only_w=True) #jx
                    dz1 = np.diff(zw1, axis=0) #jx 
                    dKdz1 = np.diff(AKs1, axis=0)/dz1
                    dKdzf1 = dKdz1[Maskr3]
                    
            z1 = ds1['zeta'][0,:,:].values
            zf0 = zf1.copy()
            zf1 = z1[Maskr]
            if verbose:
                print('   > Prepare subsequent fields for tree %0.4f sec' % (time()-tt0))
        ds0.close()
        ds1.close()

        if counter_his == 0:
            if trim_loc == True:
                # remove points on land
                xy = np.array((plon0,plat0)).T
                #print(' ++ making pmask')
                #print(maskr.shape)
                pmask = maskr[xyT_rho_un.query(xy, workers=-1)[1]]
                #print(pmask)
                # keep only points with pmask >= maskr_crit
                pcond = pmask >= maskr_crit
                #print(pcond.sum())
                plon = plon0[pcond]
                plat = plat0[pcond]
                pcs = pcs0[pcond]
            else:
                plon = plon0.copy()
                plat = plat0.copy()
                pcs = pcs0.copy()
            # create result arrays
            NP = len(plon)
            P = dict()
            for vn in plist_main:
                # NOTE: output is packed in a dict P of arrays ordered as [time, particle]
                P[vn] = np.nan * np.ones((NTS,NP))

            # write initial positions to the results arrays
            P['lon'][it0,:] = plon
            P['lat'][it0,:] = plat
            if surface == True:
                pcs[:] = S['Cs_r'][-1]
            P['cs'][it0,:] = pcs
            
            for vn in tracer_list:
                P[vn][it0,:] = get_VR(trf0_dict[vn], trf1_dict[vn], plon, plat, pcs, 0, surface)
                
            V = get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, 0, surface)
            ZH = get_zh(zf0,zf1,hf, plon, plat, 0)
            P['u'][it0,:] = V[:,0]
            P['v'][it0,:] = V[:,1]
            P['w'][it0,:] = V[:,2]
            P['zeta'][it0,:] = ZH[:,0]
            P['h'][it0,:] = ZH[:,1]
            P['z'][it0,:] = pcs * ZH.sum(axis=1) + ZH[:,0]
            
        # do the particle tracking for a single pair of history files in ndiv steps
        it1 = it0
        tt00 = time()
        for nd in range(ndiv):
            fr0 = nd/ndiv
            fr1 = (nd + 1)/ndiv
            frmid = (fr0 + fr1)/2
            
            if TR['no_advection'] == False:
                # RK4 integration
                V0 = get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, fr0, surface)
                ZH0 = get_zh(zf0,zf1,hf, plon, plat, fr0)
                if jx == False:
                    plon1, plat1, pcs1 = update_position(dxg, dyg, maskr, V0, ZH0, S, delt/2,
                                                        plon, plat, pcs, surface)
                else:
                    plon1, plat1 = update_position_lonlat(dxg, dyg, maskr, V0, delt/2, plon, plat)
                    ZH0_new = get_zh(zf0,zf1,hf, plon1, plat1, fr0)
                    pcs1 = update_position_z(V0, ZH0, ZH0_new, delt/2, pcs, surface)
                
                
                V1 = get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon1, plat1, pcs1, frmid, surface)
                ZH1 = get_zh(zf0,zf1,hf, plon1, plat1, frmid)

                if jx == False:
                    plon2, plat2, pcs2 = update_position(dxg, dyg, maskr, V1, ZH1, S, delt/2,
                                                        plon, plat, pcs, surface)
                else:
                    plon2, plat2 = update_position_lonlat(dxg, dyg, maskr, V1, delt/2, plon, plat)
                    ZH1_new = get_zh(zf0,zf1,hf, plon2, plat2, frmid)
                    pcs2 = update_position_z(V1, ZH1, ZH1_new, delt/2, pcs, surface)
                
                
                V2 = get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon2, plat2, pcs2, frmid, surface)
                ZH2 = get_zh(zf0,zf1,hf, plon2, plat2, frmid)

                if jx == False:
                    plon3, plat3, pcs3 = update_position(dxg, dyg, maskr, V2, ZH2, S, delt,
                                                        plon, plat, pcs, surface)
                else:
                    plon3, plat3 = update_position_lonlat(dxg, dyg, maskr, V2, delt, plon, plat)
                    ZH2_new = get_zh(zf0,zf1,hf, plon3, plat3, frmid)
                    pcs3 = update_position_z(V2, ZH2, ZH2_new, delt, pcs, surface)
                
                
                V3 = get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon3, plat3, pcs3, fr1, surface)
                ZH3 = get_zh(zf0,zf1,hf, plon3, plat3, fr1) # not needed?
                
                # add windage, calculated from the middle time
                if (surface == True) and (windage > 0):
                    Vwind3 = get_wind(Uwindf0, Uwindf1, Vwindf0, Vwindf1, plon, plat, frmid, windage)
                else:
                    Vwind3 = np.zeros((NP,3))
                
                # add sinking speed
                if TR['sink'] != 0:
                    Vsink3 = np.zeros((NP,3))
                    Vsink3[:,2] = -TR['sink'] / 86400
                else:
                    Vsink3 = np.zeros((NP,3))
                
                # stay close to a set depth TR['stay'] (or the bottom if it is shallower)
                if TR['stay'] != 0:
                    stay_H = ZH1.sum(axis=1)
                    stay_h = ZH1[:,1]
                    stay_z = stay_H*pcs2 + ZH1[:,0]
                    # head for the shallower of the target depth or the bottom
                    stay_target_depth = np.minimum(stay_h, TR['stay'])
                    stay_target_z = - stay_target_depth
                    stay_dz = stay_z - stay_target_z
                    # Note TR['stay'] is a depth, assumed positive down, so the
                    # signs can be confusing because z is positive up. That is why we introduce
                    # stay_target_z, which is positive up
                    Vstay3 = np.zeros((NP,3))
                    # move halfway to the desired depth each timestep
                    Vstay3[:,2] = - stay_dz / (2 * delt)
                else:
                    Vstay3 = np.zeros((NP,3))
                
                if jx == False:
                    plon, plat, pcs = update_position(dxg, dyg, maskr, (V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3 + Vsink3 + Vstay3, (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6, S, delt, plon, plat, pcs, surface)
                else:
                    plon, plat = update_position_lonlat(dxg, dyg, maskr, (V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3 + Vsink3 + Vstay3, delt, plon, plat)
                    ZH_new = get_zh(zf0,zf1,hf, plon, plat, fr1)  # use fr1? yes I think so
                    pcs = update_position_z((V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3 + Vsink3 + Vstay3,
                                            ZH0,
                                            ZH_new, delt, pcs, surface)
                
                
            elif TR['no_advection'] == True:
                V3 = np.zeros((NP,3))
                ZH3 = get_zh(zf0,zf1,hf, plon, plat, frmid)
                                              
            # add turbulence to vertical position change (advection already added above)
            if turb == True:
                # pull values of VdAKs and add up to 3-dimensions
                VdAKs = get_dAKs_new(dKdzf0, dKdzf1, plon, plat, pcs, frmid)
                VdAKs3 = np.zeros((NP,3))
                VdAKs3[:,2] = VdAKs
                # update position advecting vertically with 1/2 of AKs gradient
                ZH = get_zh(zf0,zf1,hf, plon, plat, frmid)
                #plon_junk, plat_junk, pcs_half = update_position(dxg, dyg, maskr,
                #                VdAKs3/2, ZH, S, delt/2, plon, plat, pcs, surface)
                pcs_half = update_position_z(VdAKs3/2, ZH, ZH, delt/2, pcs, surface)
                # get AKs at this height, and thence the turbulent perturbation velocity
                Vturb = get_turb(VdAKs, AKsf0, AKsf1, delt, plon, plat, pcs_half, frmid)
                Vturb3 = np.zeros((NP,3))
                Vturb3[:,2] = Vturb
                # update vertical position for real
                if jx == False:
                    plon_junk, plat_junk, pcs = update_position(dxg, dyg, maskr, Vturb3, ZH, S, delt,
                                                        plon, plat, pcs, surface)
                else:
                    pcs = update_position_z(Vturb3, ZH_new, ZH_new, delt, pcs, surface)  #jx

            ihr = nd + 1 # number of fractions 1/ndiv into the hour
            nihr = int(ndiv/TR['sph']) # number of fractions 1/ndiv between saves
            if np.mod(ihr,nihr) == 0:
                it1 += 1
                # write positions to the results arrays
                P['lon'][it1,:] = plon
                P['lat'][it1,:] = plat
                if surface == True:
                    pcs[:] = S['Cs_r'][-1]
                P['cs'][it1,:] = pcs
                
                for vn in tracer_list:
                    P[vn][it1,:] = get_VR(trf0_dict[vn], trf1_dict[vn], plon, plat, pcs, fr1, surface)

                P['u'][it1,:] = V3[:,0]
                P['v'][it1,:] = V3[:,1]
                P['w'][it1,:] = V3[:,2]
                P['zeta'][it1,:] = ZH3[:,0]
                P['h'][it1,:] = ZH3[:,1]
                P['z'][it1,:] = pcs * ZH3.sum(axis=1) + ZH3[:,0]
        if verbose:
            print('   > RK4 integration took %0.4f sec' % (time()-tt00))
        
    # and save the time vector (seconds in whatever the model reports)
    P['ot'] = rot_save

    return P

def update_position_lonlat(dxg, dyg, maskr, V, dt_sec, plon, plat):
    # find the new position
    Plon = plon.copy()
    Plat = plat.copy()
    # This next step is the actual particle displacement, just done
    # by velocity times a time interval, giving displacements in meters.
    # Each row is a different particle and the columns are (x,y,z) = [0,1,2]
    dX_m = V*dt_sec
    # Horizontal advection
    per_m = zfun.earth_rad(Plat)
    clat = np.cos(np.pi*Plat/180.)
    pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
    pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
    Plon += pdx_deg
    Plat += pdy_deg
    
    # Keep particles from being trapped on land
    # ## * minimum grid sizes used when particles approach land boundaries
    # Experiments with "trap0" to explore trapping in the Skokomish
    # showed that ## = 0.5 is a reasonable choice.
    xy = np.array((Plon,Plat)).T
    pmask = maskr[xyT_rho_un.query(xy, workers=-1)[1]]
    pcond = pmask < maskr_crit # a Boolean mask
    if len(pcond) > 0:
        # these randint calls give random vectors of -1,0,1 (note the 2!)
        rix = np.random.randint(-1,2,len(plon))
        riy = np.random.randint(-1,2,len(plon))
        Plon[pcond] = plon[pcond] + 0.5*rix[pcond]*dxg
        Plat[pcond] = plat[pcond] + 0.5*riy[pcond]*dyg
        
    # move any particles on land to the middle of the nearest good rho point.
    xy = np.array((Plon,Plat)).T
    pmask = maskr[xyT_rho_un.query(xy, workers=-1)[1]]
    pcond = pmask < maskr_crit # a Boolean mask
    if len(pcond) > 0:
        Plon_NEW = lonrf[xyT_rho.query(xy, workers=-1)[1]]
        Plat_NEW = latrf[xyT_rho.query(xy, workers=-1)[1]]
        Plon[pcond] = Plon_NEW[pcond]
        Plat[pcond] = Plat_NEW[pcond]
    return Plon, Plat

def update_position_z(V, ZH0, ZH1, dt_sec, pcs0, surface):
    # find the new position in the vertical
    Pcs = pcs0.copy()
    # vertical particle displacement
    dX_m = V*dt_sec
    # Vertical advection
    H0 = ZH0.sum(axis=1); # ZH0 is the current total depth, ZH1 is the total depth after updating plon and plat
    # NOTE: The first column of ZH is the
    # zeta of all particles, and the second is bottom depth h (a positive number),
    # so the H we calculate above is the total water column thickness.
    # The particle z-positions are given by H*pcs + ZH[:,0]
    z_par0 = H0*pcs0 + ZH0[:,0]
    z_par1 = z_par0 + dX_m[:,2]
    zeta1 = [x[0] for x in ZH1]
    pcs1 = (z_par1 - zeta1)/ ZH1.sum(axis=1)
 #   pdz_s = dX_m[:,2]/H
    # Reflective upper and lower boundary conditions
    if surface == False:
    #    Pcs_orig = Pcs.copy()
    #    Pcs += pdz_s
        Pcs = pcs1.copy()
        # Check for bad pcs values.  This was happening a lot becasue of
        # some bugs on the turbulence code that would return bad vertical velocities
        # but since those were fixed this "bad_pcs" diagnostic is all zeros.
        pcs_mask = np.isnan(Pcs)
        if sum(pcs_mask) > 0:
            print('Warning: bad pcs values')
            Pcs[pcs_mask] = Pcs_orig[pcs_mask]
        # Enforce limits on cs using reflection.  We use the remainder function to
        # account for cases where the vertical advection may have moved
        # particles more than 2*H
        hit_top = Pcs > 0
        Pcs[hit_top] = - np.remainder(Pcs[hit_top],1)
        hit_bottom = Pcs < -1
        Pcs[hit_bottom] = -1 - np.remainder(Pcs[hit_bottom],-1)
        # and finally enforce more limits if needed
        Pcs[Pcs < -1] = -1
        Pcs[Pcs > 0] = 0
    else:
        Pcs[:] = 0 # surface trapped

    return Pcs
    
def update_position(dxg, dyg, maskr, V, ZH, S, dt_sec, plon, plat, pcs, surface):
    
    # find the new position
    Plon = plon.copy()
    Plat = plat.copy()
    Pcs = pcs.copy()
    # This next step is the actual particle displacement, just done
    # by velocity times a time interval, giving displacements in meters.
    # Each row is a different particle and the columns are (x,y,z) = [0,1,2]
    dX_m = V*dt_sec
    # Horizontal advection
    per_m = zfun.earth_rad(Plat)
    clat = np.cos(np.pi*Plat/180.)
    pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
    pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
    Plon += pdx_deg
    Plat += pdy_deg
    
    # Keep particles from being trapped on land
    # ## * minimum grid sizes used when particles approach land boundaries
    # Experiments with "trap0" to explore trapping in the Skokomish
    # showed that ## = 0.5 is a reasonable choice.
    xy = np.array((Plon,Plat)).T
    pmask = maskr[xyT_rho_un.query(xy, workers=-1)[1]]
    pcond = pmask < maskr_crit # a Boolean mask
    if len(pcond) > 0:
        # these randint calls give random vectors of -1,0,1 (note the 2!)
        rix = np.random.randint(-1,2,len(plon))
        riy = np.random.randint(-1,2,len(plon))
        Plon[pcond] = plon[pcond] + 0.5*rix[pcond]*dxg
        Plat[pcond] = plat[pcond] + 0.5*riy[pcond]*dyg
        
    # move any particles on land to the middle of the nearest good rho point.
    xy = np.array((Plon,Plat)).T
    pmask = maskr[xyT_rho_un.query(xy, workers=-1)[1]]
    pcond = pmask < maskr_crit # a Boolean mask
    if len(pcond) > 0:
        Plon_NEW = lonrf[xyT_rho.query(xy, workers=-1)[1]]
        Plat_NEW = latrf[xyT_rho.query(xy, workers=-1)[1]]
        Plon[pcond] = Plon_NEW[pcond]
        Plat[pcond] = Plat_NEW[pcond]
        
    # Start of Vertical advection
    H = ZH.sum(axis=1)
    # NOTE: The first column of ZH is the
    # zeta of all particles, and the second is bottom depth h (a positive number),
    # so the H we calculate above is the total water column thickness.
    # The particle z-positions are given by H*pcs + ZH[:,0]
    pdz_s = dX_m[:,2]/H
    # Reflective upper and lower boundary conditions
    if surface == False:
        Pcs_orig = Pcs.copy()
        Pcs += pdz_s
        # Check for bad pcs values.  This was happening a lot becasue of
        # some bugs on the turbulence code that would return bad vertical velocities
        # but since those were fixed this "bad_pcs" diagnostic is all zeros.
        pcs_mask = np.isnan(Pcs)
        if sum(pcs_mask) > 0:
            print('Warning: bad pcs values')
            Pcs[pcs_mask] = Pcs_orig[pcs_mask]
        # Enforce limits on cs using reflection.  We use the remainder function to
        # account for cases where the vertical advection may have moved
        # particles more than 2*H 
        hit_top = Pcs > 0
        Pcs[hit_top] = - np.remainder(Pcs[hit_top],1)
        hit_bottom = Pcs < -1
        Pcs[hit_bottom] = -1 - np.remainder(Pcs[hit_bottom],-1)
        # and finally enforce more limits if needed
        Pcs[Pcs < -1] = -1
        Pcs[Pcs > 0] = 0
    else:
        Pcs[:] = 0 # surface trapped

    return Plon, Plat, Pcs

def get_vel(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, frac, surface):
    # Get the velocity at all points, at an arbitrary time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1, 0 <= frac <= 1.
    # NOTE: with ndiv=1 this gets called 4 times per hour, or 96 times per day.
    NP = len(plon)
    V = np.zeros((NP,3))
    if surface == True:
        xy = np.array((plon,plat)).T
        # use workers=-1 to use all available cores
        ui0 = uf0[xyT_u.query(xy, workers=-1)[1]]
        vi0 = vf0[xyT_v.query(xy, workers=-1)[1]]
        ui1 = uf1[xyT_u.query(xy, workers=-1)[1]]
        vi1 = vf1[xyT_v.query(xy, workers=-1)[1]]
        ui = (1 - frac)*ui0 + frac*ui1
        vi = (1 - frac)*vi0 + frac*vi1
        V[:,0] = ui
        V[:,1] = vi
    else:
        xys = np.array((plon,plat,pcs)).T
        ui0 = uf0[xyzT_u.query(xys, workers=-1)[1]]
        vi0 = vf0[xyzT_v.query(xys, workers=-1)[1]]
        wi0 = wf0[xyzT_w.query(xys, workers=-1)[1]]
        ui1 = uf1[xyzT_u.query(xys, workers=-1)[1]]
        vi1 = vf1[xyzT_v.query(xys, workers=-1)[1]]
        wi1 = wf1[xyzT_w.query(xys, workers=-1)[1]]
        ui = (1 - frac)*ui0 + frac*ui1
        vi = (1 - frac)*vi0 + frac*vi1
        wi = (1 - frac)*wi0 + frac*wi1
        V[:,0] = ui
        V[:,1] = vi
        V[:,2] = wi
    # The line below is to account for wet-dry cells.
    V[np.isnan(V)] = 0.0
    return V
    
def get_zh(zf0,zf1,hf, plon, plat, frac):
    # Get zeta and h at all points, at an arbitrary time between two saves
    NP = len(plon)
    xy = np.array((plon,plat)).T
    zi0 = zf0[xyT_rho.query(xy, workers=-1)[1]]
    zi1 = zf1[xyT_rho.query(xy, workers=-1)[1]]
    hi = hf[xyT_rho.query(xy, workers=-1)[1]]
    zi = (1 - frac)*zi0 + frac*zi1
    ZH = np.zeros((NP,2))
    ZH[:,0] = zi
    ZH[:,1] = hi
    return ZH
    
def get_VR(tf0,tf1, plon, plat, pcs, frac, surface):
    # Get a variable on the z_rho grid at all points.
    if surface == True:
        xy = np.array((plon,plat)).T
        ti0 = tf0[xyT_rho.query(xy, workers=-1)[1]]
        ti1 = tf1[xyT_rho.query(xy, workers=-1)[1]]
    else:
        xys = np.array((plon,plat,pcs)).T
        ti0 = tf0[xyzT_rho.query(xys, workers=-1)[1]]
        ti1 = tf1[xyzT_rho.query(xys, workers=-1)[1]]
    ti = (1 - frac)*ti0 + frac*ti1
    return ti
    
def get_wind(Uwindf0, Uwindf1, Vwindf0, Vwindf1, plon, plat, frac, windage):
    # creates the windage correction to the surface velocity (u,v only)
    NP = len(plon)
    Vwind3 = np.zeros((NP,3))
    xy = np.array((plon,plat)).T
    Uwind00 = Uwindf0[xyT_rho.query(xy, workers=-1)[1]]
    Uwind11 = Uwindf1[xyT_rho.query(xy, workers=-1)[1]]
    Uwind = (1 - frac)*Uwind00 + frac*Uwind11
    Vwind00 = Vwindf0[xyT_rho.query(xy, workers=-1)[1]]
    Vwind11 = Vwindf1[xyT_rho.query(xy, workers=-1)[1]]
    Vwind = (1 - frac)*Vwind00 + frac*Vwind11
    Vwind3[:,0] = windage*Uwind
    Vwind3[:,1] = windage*Vwind
    return Vwind3
    
def get_AKs(AKsf, plon, plat, pcs):
    # Get AKs at all points, at one time.
    xys = np.array((plon,plat,pcs)).T
    AKsi = AKsf[xyzT_w.query(xys, workers=-1)[1]]
    return AKsi
    
def get_dAKs_new(dKdzf0, dKdzf1, plon, plat, pcs, frac):
    xys = np.array((plon,plat,pcs)).T
    dKdzi0 = dKdzf0[xyzT_rho.query(xys, workers=-1)[1]]
    dKdzi1 = dKdzf1[xyzT_rho.query(xys, workers=-1)[1]]
    dKdzi = (1 - frac)*dKdzi0 + frac*dKdzi1
    return dKdzi

def get_turb(dAKs, AKsf0, AKsf1, delta_t, plon, plat, pcs, frac):
    # get the vertical turbulence correction components
    V0 = get_AKs(AKsf0, plon, plat, pcs)
    V1 = get_AKs(AKsf1, plon, plat, pcs)
    # create weighted average diffusivity
    Vave = (1 - frac)*V0 + frac*V1
    # turbulence calculation from Banas, MacCready, and Hickey (2009)
    # w_turbulence = rand*sqrt(2K/dt) + dK/dz
    # rand = random array with normal distribution
    rand = np.random.standard_normal(len(V0))
    V = rand*np.sqrt(2*Vave/delta_t) + dAKs
    return V

def get_fn_list(idt, Ldir):
    # Gets history files for 1 day only:
    # ocean_his_0025.nc from the day before through ocean_his_0025.nc of this day.
    ds0 = idt.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds0)
    return fn_list
    
