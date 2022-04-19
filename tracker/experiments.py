"""
This is where you set the run initial condition using get_ic()
based on an experiment name passed by the calling code.

Thre are also some utility functions useful for making different
common release patterns.

"""

import numpy as np
    
def get_ic(TR):
    # routines to set particle initial locations, all numpy arrays
    
    # NOTE: "pcs" refers to fractional depth, and goes linearly from -1 to 0
    # between the local bottom and free surface.  It is how we keep track of
    # vertical position, only converting to z-position when needed.
    
    exp_name = TR['exp_name']
    gridname = TR['gridname']
    fn00 = TR['fn00']
        
    if exp_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
    
    if exp_name == 'ai0': # Mid-Admiralty Inlet
        lonvec = np.array([-122.6])
        latvec = np.array([48])
        # These are: (Slope off JdF, Middle of JdF, Whidbey Basin)
        pcs_vec = np.linspace(-1,-0.9,num=1000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'vmix': # three vertical profiles to test mixing
        # use with the new flag: -no_advection True, so a full command would be
        # python tracker.py -exp vmix -3d True -clb True -no_advection True
        lonvec = np.array([-125.35, -124.0, -122.581])
        latvec = np.array([47.847, 48.3, 48.244])
        # These are: (Slope off JdF, Middle of JdF, Whidbey Basin)
        pcs_vec = np.linspace(-1,0,num=4000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'dmMerhab':
        nyp = 7
        x0 = -126; x1 = -125; y0 = 48; y1 = 49
        clat_1 = np.cos(np.pi*(np.mean([y0, y1]))/180)
        xyRatio = clat_1 * (x1 - x0) / (y1 - y0)
        lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        latvec = np.linspace(y0, y1, nyp)
        lonmat_1, latmat_1 = np.meshgrid(lonvec, latvec)
        #
        x0 = -125.2; x1 = -124.2; y0 = 44; y1 = 45
        clat_2 = np.cos(np.pi*(np.mean([y0, y1]))/180)
        xyRatio = clat_2 * (x1 - x0) / (y1 - y0)
        lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        latvec = np.linspace(y0, y1, nyp)
        lonmat_2, latmat_2 = np.meshgrid(lonvec, latvec)
        lonmat = np.concatenate((lonmat_1.flatten(), lonmat_2.flatten()))
        latmat = np.concatenate((latmat_1.flatten(), latmat_2.flatten()))
        #
        plon00 = lonmat.flatten(); plat00 = latmat.flatten()
        pcs00 = np.zeros(plon00.shape)
        
    elif exp_name == 'full': # the whole domain of cas6, with some edges trimmed
        # used by drifters0
        lonvec = np.linspace(-129, -122, 30)
        latvec = np.linspace(43, 51, 60)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'PS': # nominally Puget Sound
        # used by drifters0
        lonvec = np.linspace(-123.6, -122, 50)
        latvec = np.linspace(47, 49, 100)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    return plon00, plat00, pcs00
    
def ic_from_meshgrid(lonvec, latvec, pcs_vec):
    # First create three vectors of initial locations (as done in some cases above).
    # plat00 and plon00 should be the same length, and the length of pcs00 is
    # as many vertical positions you have at each lat, lon
    # (expressed as fraction of depth -1 < pcs < 0).
    # Then we create full output vectors (each has one value per point).
    # This code takes each lat, lon location and then assigns it to NSP points
    # corresponding to the vector of pcs values.
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon vectors')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00
    
def ic_from_list(lonvec, latvec, pcs_vec):
    # Like ic_from_meshgrid() but treats the lon, lat lists like lists of mooring locations.
    plon_vec = lonvec
    plat_vec = latvec
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon lists')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00
    
def ic_from_TEFsegs(fn00, gridname, seg_list, DZ, NPmax=10000):
    import pickle
    import sys
    # select the indir
    from lo_tools import Lfun, zrfun
    Ldir = Lfun.Lstart()
    indir = Ldir['LOo'] / 'tef' / ('volumes_' + gridname)
    # load data
    j_dict = pickle.load(open(indir / 'j_dict.p', 'rb'))
    i_dict = pickle.load(open(indir / 'i_dict.p', 'rb'))
    G = zrfun.get_basic_info(fn00, only_G=True)
    h = G['h']
    xp = G['lon_rho']
    yp = G['lat_rho']
    plon_vec = np.array([])
    plat_vec = np.array([])
    hh_vec = np.array([])
    for seg_name in seg_list:
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        # untested 2021.10.05
        hh_vec = np.append(hh_vec, h[jjj,iii])
        plon_vec = np.append(plon_vec, xp[jjj,iii])
        plat_vec = np.append(plat_vec, yp[jjj,iii])
        # ji_seg = ji_dict[seg_name]
        # for ji in ji_seg:
        #     plon_vec = np.append(plon_vec, xp[ji])
        #     plat_vec = np.append(plat_vec, yp[ji])
        #     hh_vec = np.append(hh_vec, h[ji])
    plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
    for ii in range(len(plon_vec)):
        x = plon_vec[ii]
        y = plat_vec[ii]
        hdz = DZ*np.floor(hh_vec[ii]/DZ) # depth to closest DZ m (above the bottom)
        if hdz >= DZ:
            zvec = np.arange(-hdz,DZ,DZ) # a vector that goes from -hdz to 0 in steps of DZ m
            svec = zvec/hh_vec[ii]
            ns = len(svec)
            if ns > 0:
                plon00 = np.append(plon00, x*np.ones(ns))
                plat00 = np.append(plat00, y*np.ones(ns))
                pcs00 = np.append(pcs00, svec)
    # subsample the I.C. vectors to around max length around NPmax
    NP = len(plon00)
    print(len(plon00))
    nstep = max(1,int(NP/NPmax))
    plon00 = plon00[::nstep]
    plat00 = plat00[::nstep]
    pcs00 = pcs00[::nstep]
    print(len(plon00))
    sys.stdout.flush()
    return plon00, plat00, pcs00
    