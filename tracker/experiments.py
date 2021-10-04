"""
This is where you set the run "gtagex" and the initial condition
based on an experiment name passed by the calling code.

"""

import numpy as np
import sys

def get_exp_info(exp_name):
    
    # Defaults
    # you could override these using "if exp_name == ..."
    gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
    
    if exp_name == 'EJdF3d':
        ex_name = 'lo8da'
    
    EI = {}
    EI['exp_name'] = exp_name # tracker experiment name
    # ROMS names
    EI['gridname'] = gridname
    EI['tag'] = tag
    EI['ex_name'] = ex_name
    EI['gtagex'] = gridname + '_' + tag + '_' + ex_name
    return EI
    
def get_ic(EI, fn00):
    # routines to set particle initial locations, all numpy arrays
    
    # NOTE: "pcs" refers to fractional depth, and goes linearly from -1 to 0
    # between the local bottom and free surface.  It is how we keep track of
    # vertical position, only converting to z-position when needed.
    
    exp_name = EI['exp_name']
    
    if exp_name == 'p5_PS':
        rloc_dict = {'Whidbey NW': (-122.749, 48.2927),
                    'Penn Cove': (-122.695, 48.2309)}
        rnp = 100
        count = 0
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]; ry = rloc_dict[rloc][1]
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy(); plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        pcs00 = np.zeros_like(plon00)
    
    elif exp_name == 'dcrab1':
        rloc_dict = {'South Sound': (-122.7, 47.11),
                    'Hood Canal': (-122.86, 47.67),
                    'Whidbey': (-122.6, 48.25)}
        rnp = 100
        count = 0
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]; ry = rloc_dict[rloc][1]
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy(); plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        pcs00 = np.zeros_like(plon00)
        
    elif exp_name == 'dcrab2':
        rloc_dict = {'Grays Harbor': (-124.3, 46.9)}
        rnp = 100
        count = 0
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]; ry = rloc_dict[rloc][1]
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy(); plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        pcs00 = np.zeros_like(plon00)
    
    elif exp_name == 'hat':
        rloc_dict = {'Near Hat Island': (-122.27, 48)}
        rnp = 100
        count = 0
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]; ry = rloc_dict[rloc][1]
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.006*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy(); plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        pcs00 = np.zeros_like(plon00)
    
    elif exp_name == 'full': # the whole domain of cas6, with some edges trimmed
        lonvec = np.linspace(-129, -122, 30)
        latvec = np.linspace(43, 51, 60)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'PS': # nominally Puget Sound
            lonvec = np.linspace(-123.6, -122, 50)
            latvec = np.linspace(47, 49, 100)
            pcs_vec = np.array([0])
            plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
                
    elif exp_name == 'p5_merhab':
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
        pcs00 = np.zeros_like(plon00)
        
    elif exp_name == 'HC3d': # Hood Canal 3d using rectangle
        # This fills a volume defined by a rectangular lon, lat region
        # with particles spaced every DZ m from the surface to near the bottom
        lonvec = np.linspace(-123.15, -122.9, 30)
        latvec = np.linspace(47.45, 47.65, 30)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        import zfun
        import zrfun
        G, S, T = zrfun.get_basic_info(fn00)
        hh = zfun.interp2(lonmat, latmat, G['lon_rho'], G['lat_rho'],G['h'])
        plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        hh_vec = hh.flatten()
        DZ = 10
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
                    
    elif exp_name == 'InnerHC3d': # Hood Canal 3d using TEF segments
        seg_list = ['H3','H4','H5','H6','H7','H8']
        DZ = 5
        plon00, plat00, pcs00 = ic_from_TEFsegs(fn00, EI['gridname'], seg_list, DZ)

    elif exp_name == 'AllHC3d': # Hood Canal 3d using TEF segments
        seg_list = ['H'+str(item) for item in range(1,9)]
        DZ = 2
        plon00, plat00, pcs00 = ic_from_TEFsegs(fn00, EI['gridname'], seg_list, DZ, NPmax=100000)
        
    elif exp_name == 'EJdF3d': # eastern JdF 3d using a TEF segment
        seg_list = ['J4']
        DZ = 5
        plon00, plat00, pcs00 = ic_from_TEFsegs(fn00, EI['gridname'], seg_list, DZ, NPmax=10000)

    elif exp_name == 'tn0': # Tacoma Narrows region
        #to test headland jumping
        lonvec = np.linspace(-122.65, -122.45, 30)
        latvec = np.linspace(47.2, 47.35, 30)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'eddy0': # Same as particulator JdF eddy release
        lonvec = np.linspace(-125.6, -125.2, 20)
        latvec = np.linspace(48.4, 48.6, 20)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'skok': # head of the Skokomish River
        # to test trapping in rivers
        lonvec = np.linspace(-123.171, -123.163, 10)
        latvec = np.linspace(47.306, 47.312, 10)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'vmix': # three vertical profiles to test mixing
        # use with the new flag: -no_advection True, so a full command would be
        # python tracker.py -exp vmix -3d True -clb True -no_advection True
        lonvec = np.array([-125.35, -124.0, -122.581])
        latvec = np.array([47.847, 48.3, 48.244])
        # These are: (Slope off JdF, Middle of JdF, Whidbey Basin)
        pcs_vec = np.linspace(-1,0,num=4000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'girton1': # single water column for James Girton
        # 1) Date: 6 Nov 2017, 20:47 GMT
        # Location: 47 42.812’ N, 122 25.055’ W
        # Duration: 2.3 hours
        # 2) Date: 2 Sep 2020, 18:00 GMT
        # Location: SAME
        # Duration: 4 hours
        lonvec = np.array([-(122 + 25.055/60)])
        latvec = np.array([47 + 42.812/60])
        pcs_vec = np.linspace(-1,0,num=4000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'leathermanWP': # single water column for Marissa Leatherman
        # approximates West Point Sewage Treatment Plant outfall
        lonvec = np.array([-122.46])
        latvec = np.array([47.662])
        pcs_vec = np.linspace(-.6,-.4,num=1000) # all around mid-depth
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'leathermanNavy': # single water column for Marissa Leatherman
        # approximates Puget Sound Naval Shipyard
        lonvec = np.array([-122.64])
        latvec = np.array([47.55])
        pcs_vec = np.zeros(1000) # surface
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    elif exp_name == 'leathermanPuy': # single water column for Marissa Leatherman
        # approximates Puyallup River
        lonvec = np.array([-122.43])
        latvec = np.array([47.27])
        pcs_vec = np.linspace(-1,0,num=1000) # all depths
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
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
    # select the indir
    import Lfun
    Ldir = Lfun.Lstart()
    indir = Ldir['LOo'] + 'tef2/volumes_' + gridname + '/'
    # load data
    ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))
    import zrfun
    G = zrfun.get_basic_info(fn00, only_G=True)
    h = G['h']
    xp = G['lon_rho']
    yp = G['lat_rho']
    plon_vec = np.array([])
    plat_vec = np.array([])
    hh_vec = np.array([])
    for seg_name in seg_list:
        ji_seg = ji_dict[seg_name]
        for ji in ji_seg:
            plon_vec = np.append(plon_vec, xp[ji])
            plat_vec = np.append(plat_vec, yp[ji])
            hh_vec = np.append(hh_vec, h[ji])
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
    