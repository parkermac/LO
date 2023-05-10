# -*- coding: utf-8 -*-
"""
Module of functions specific to ROMS.

2021.08.05 Recoded using xarray instead of netCDF4.
"""
import xarray as xr
import numpy as np
import pandas as pd
from lo_tools import Lfun
import sys
import pickle

def get_basic_info(fn, only_G=False, only_S=False, only_T=False):
    """
    Gets grid, vertical coordinate, and time info from a ROMS NetCDF
    history file with full name 'fn'
    Input: the filename (with path if needed)
    Output: dicts G, S, and T
    Example calls:
    G, S, T = zfun.get_basic_info(fn)
    T = zfun.get_basic_info(fn, only_T=True)
    """
    ds = xr.open_dataset(fn)
    def make_G(ds):
        # get grid and bathymetry info
        g_varlist = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v',
        'lon_psi', 'lat_psi', 'mask_rho', 'mask_u', 'mask_v', 'pm', 'pn',]
        G = dict()
        for vv in g_varlist:
            G[vv] = ds[vv].values
        G['DX'] = 1/G['pm']
        G['DY'] = 1/G['pn']
        G['M'], G['L'] = np.shape(G['lon_rho']) # M = rows, L = columns
        return G
    def make_S(ds):
        # get vertical sigma-coordinate info (vectors are bottom to top)
        s_varlist = ['s_rho', 's_w', 'hc', 'Cs_r', 'Cs_w', 'Vtransform']
        S = dict()
        for vv in s_varlist:
            S[vv] = ds[vv].values
        S['N'] = len(S['s_rho']) # number of vertical levels
        return S
    def make_T(ds):
        # returns two single values, one a datatime, and one a float
        ot = ds.ocean_time.values # an array with dtype='datetime64[ns]'
        dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
        dt = dti.to_pydatetime() # an array of datetimes
        T = dict()
        T['dt'] = dt[0] # a datetime object
        T['ocean_time'] = Lfun.datetime_to_modtime(dt[0]) # a float, "seconds since..."
        return T
    # return results
    if only_G:
        return make_G(ds)
    elif only_S:
        return make_S(ds)
    elif only_T:
        return make_T(ds)
    else:
        return make_G(ds), make_S(ds), make_T(ds)
    ds.close()

def get_z(h, zeta, S, only_rho=False, only_w=False):
    """
    Used to calculate the z position of fields in a ROMS history file

    Input: arrays h (bathymetry depth) and zeta (sea surface height)
    which must be the same size, and dict S created by get_basic_info()

    Output: 3-D arrays of z_rho and z_w

    NOTE: one foible is that if you input arrays of h and zeta that are
    vectors of length VL, the output array (e.g. z_rho) will have size (N, VL)
    (i.e. it will never return an array with size (N, VL, 1), even if (VL, 1) was
    the input shape).  This is a result of the initial and final squeeze calls.
    """
    # input error checking
    if ( (not isinstance(h, np.ndarray))
        or (not isinstance(zeta, (np.ndarray, np.ma.core.MaskedArray))) ):
        print('WARNING from get_z(): Inputs must be numpy arrays')
    if not isinstance(S, dict):
        print('WARNING from get_z(): S must be a dict')
    # number of vertical levels
    N = S['N']
    # remove singleton dimensions
    h = h.squeeze()
    zeta = zeta.squeeze()
    # ensure that we have enough dimensions
    h = np.atleast_2d(h)
    zeta = np.atleast_2d(zeta)
    # check that the dimensions are the same
    if h.shape != zeta.shape:
        print('WARNING from get_z(): h and zeta must be the same shape')
    M, L = h.shape
    def make_z_rho(h, zeta, S, N, M, L):
        # rho
        # create some useful arrays
        csr = S['Cs_r']
        csrr = csr.reshape(N, 1, 1).copy()
        Cs_r = np.tile(csrr, [1, M, L])
        H_r = np.tile(h.reshape(1, M, L).copy(), [N, 1, 1])
        Zeta_r = np.tile(zeta.reshape(1, M, L).copy(), [N, 1, 1])
        if S['hc'] == 0: # if hc = 0 the transform is simpler (and faster)
            z_rho = H_r*Cs_r + Zeta_r + Zeta_r*Cs_r
        elif S['hc'] != 0: # need to calculate a few more useful arrays
            sr = S['s_rho'] # PM edit 2019.01.24
            srr = sr.reshape(N, 1, 1).copy()
            S_rho = np.tile(srr, [1, M, L])
            Hc_r = np.tile(S['hc'], [N, M, L])
            if S['Vtransform'] == 1:
                zr0 = (S_rho - Cs_r) * Hc_r + Cs_r*H_r
                z_rho = zr0 + Zeta_r * (1 + zr0/H_r)
            elif S['Vtransform'] == 2:
                zr0 = (S_rho*Hc_r + Cs_r*H_r) / (Hc_r + H_r)
                z_rho = Zeta_r + (Zeta_r + H_r)*zr0
        z_rho = z_rho.squeeze()
        return z_rho
    def make_z_w(h, zeta, S, N, M, L):
        # w
        # create some useful arrays
        csw = S['Cs_w']
        csww = csw.reshape(N+1, 1, 1).copy()
        Cs_w = np.tile(csww, [1, M, L])
        H_w = np.tile(h.reshape(1, M, L).copy(), [N+1, 1, 1])
        Zeta_w = np.tile(zeta.reshape(1, M, L).copy(), [N+1, 1, 1])
        if S['hc'] == 0: # if hc = 0 the transform is simpler (and faster)
            z_w = H_w*Cs_w + Zeta_w + Zeta_w*Cs_w
        elif S['hc'] != 0: # need to calculate a few more useful arrays
            #sw = S['s_w']
            sw = S['s_w'] # PM edit 2019.01.24
            sww = sw.reshape(N+1, 1, 1).copy()
            S_w = np.tile(sww, [1, M, L])    #
            Hc_w = np.tile(S['hc'], [N+1, M, L])
            if S['Vtransform'] == 1:
                zw0 = (S_w - Cs_w) * Hc_w + Cs_w*H_w
                z_w = zw0 + Zeta_w * (1 + zw0/H_w)
            elif S['Vtransform'] == 2:
                zw0 = (S_w*Hc_w  + Cs_w*H_w) / (Hc_w + H_w)
                z_w = Zeta_w + (Zeta_w + H_w)*zw0
        z_w = z_w.squeeze()
        return z_w
    # return results
    if only_rho:
        return make_z_rho(h, zeta, S, N, M, L)
    elif only_w:
        return make_z_w(h, zeta, S, N, M, L)
    else :
        return make_z_rho(h, zeta, S, N, M, L), make_z_w(h, zeta, S, N, M, L)
    
def get_S(S_info_dict):
    """
    Code to calculate S-coordinate vectors from the parameters
    in S_COORDINATE_INFO.csv.
    Need to check this carefully against the matlab version.
    # recoded for python on 7/7/2016 from:
    # Z_scoord.m  5/21/2007  Parker MacCready
    # this creates the structure S, which would be used for example by
    # Z_s2z.m, given basic grid parameters
    # edited by DAS to include more things in S stucture
    # edited by SNG March 2011 to include all of the current available ROMS
    # stretching functions, 1-4 see:
    # https://www.myroms.org/wiki/index.php/Vertical_S-coordinate#Vertical_Stretching_Functions
    
    NOTES 2019.09.11
    (1) I checked that Cs_r and _w made by this program are identical to those which are
    given in the ROMS history files.  They are.
    (2) I also made some inquiries on the ROMS forum to make sure that the parameter 'hc' is
    being done correctly.  The short answer is that yes it is.  With Vtransform = 2 (my
    new default) it is given by Tcline from the .in file.  In older runs with Vtransform = 1
    is it min(hmin, Tcline) and this REQUIRES that Tcline be less than hmin.  Since all those
    older runs used Tcline = 0 then hc = 0.
    
    """
    S = dict()
    for item in S_info_dict.keys():
        if item in ['N', 'VSTRETCHING', 'VTRANSFORM']:
            S[item.title()] = int(S_info_dict[item])
        elif item in ['TCLINE', 'THETA_S', 'THETA_B']:
            S[item.lower()] = float(S_info_dict[item])
        else:
            pass
    N = S['N']
    Vstretching = S['Vstretching']
    Vtransform = S['Vtransform']
    tcline = S['tcline']
    theta_s = S['theta_s']
    theta_b = S['theta_b']
    hmin = 3 # a placeholder, used only for Vtransform = 1.
    if Vtransform == 1:
        hc = min(hmin,tcline)
    elif Vtransform == 2:
        hc = tcline
    S['hc'] = hc
    s_rho = (np.linspace(-(N-1), 0, N) - 0.5)/N
    s_w = np.linspace(-N, 0, N+1)/N
    S['s_rho'] = s_rho
    S['s_w'] = s_w
    if Vstretching == 1:
        if theta_s != 0:
            cff1 = 1/np.sinh(theta_s)
            cff2 = 0.5/np.tanh(0.5*theta_s)
            Cs_r = ( (1-theta_b)*cff1*np.sinh(theta_s*s_rho)
                    + theta_b*( cff2*np.tanh(theta_s*(s_rho + 0.5)) - 0.5 ) )
            Cs_w = ( (1-theta_b)*cff1*np.sinh(theta_s*s_w)
                    + theta_b*( cff2*np.tanh(theta_s*(s_w + 0.5)) - 0.5 ) )
        else:
            Cs_r = s_rho
            Cs_w = s_w
    elif Vstretching == 2:
        alpha = 1
        beta = 1
        if theta_s!=0 and theta_b!=0:
            Csur = (1-np.cosh(theta_s*s_rho))/(np.cosh(theta_s)-1)
            Cbot = ((np.sinh(theta_b*(s_rho+1)))/(np.sinh(theta_b)))-1
            u = ((s_rho+1)**alpha)*(1+(alpha/beta)*(1-((s_rho+1)**beta)))
            Cs_r = u*Csur+(1-u)*Cbot
            Csur_w = (1-np.cosh(theta_s*s_w))/(np.cosh(theta_s)-1)
            Cbot_w = ((np.sinh(theta_b*(s_w+1)))/(np.sinh(theta_b)))-1
            u_w = ((s_w+1)**alpha)*(1+(alpha/beta)*(1-((s_w+1)**beta)))
            Cs_w = u_w*Csur_w+(1-u_w)*Cbot_w
        else:
            Cs_r = s_rho
            Cs_w = s_w
    elif Vstretching == 3:
        # Geyer function for high bbl resolution in shallow applications
        gamma = 3
        Csur = -(np.log(np.cosh(gamma*abs(s_rho)**theta_s)))/np.log(np.cosh(gamma))
        Cbot = ((np.log(np.cosh(gamma*(s_rho+1)**theta_b)))/np.log(np.cosh(gamma)))-1
        mu = 0.5*(1-np.tanh(gamma*(s_rho+0.5)))
        Cs_r = mu*Cbot+(1-mu)*Csur
        Csur_w = -(np.log(np.cosh(gamma*abs(s_w)**theta_s)))/np.log(np.cosh(gamma))
        Cbot_w = ((np.log(np.cosh(gamma*(s_w+1)**theta_b)))/np.log(np.cosh(gamma)))-1
        mu_w = 0.5*(1-np.tanh(gamma*(s_w+0.5)))
        Cs_w = mu_w*Cbot_w+(1-mu_w)*Csur_w
    elif Vstretching == 4:
        # newest ROMS default as of March 2011 (theta_s between 0 and 10,
        # theta_b between 0 and 4)
        if theta_s>0:
            Cs_r = (1-np.cosh(theta_s*s_rho))/(np.cosh(theta_s)-1)
            Cs_w = (1-np.cosh(theta_s*s_w))/(np.cosh(theta_s)-1)
        elif theta_s<=0:
            Cs_r = -(s_rho**2)
            Cs_w = -(s_w**2)
        if theta_b > 0:
            Cs_r = (np.exp(theta_b*Cs_r)-1)/(1-np.exp(-theta_b))
            Cs_w = (np.exp(theta_b*Cs_w)-1)/(1-np.exp(-theta_b))
    S['Cs_r'] = Cs_r
    S['Cs_w'] = Cs_w
    return S

def make_varinfo_list():
    """
    This method pre-parses varinfo.yaml into a list for faster use by get_varinfo()
    
    Configured to use the new ROMS varinfo.yaml (1/2022), but you need a few edits
    so we use our modified version in LO_roms_source/npzd_banas.
    """
    import yaml
        
    # specify which varinfo.yaml to use
    Ldir = Lfun.Lstart()
    fn = Ldir['parent'] / 'LO_roms_source_alt' / 'varinfo' / 'varinfo.yaml'
    
    # parse into a list of dicts
    with open(fn,'r') as f:
        yaml_dict = yaml.safe_load(f)
    short_list = yaml_dict['metadata'] # a list of dicts, one per item
    
    # remove some things from the list
    short_list = [item for item in short_list if 'adjoint' not in item['field']]
    short_list = [item for item in short_list if 'tangent' not in item['field']]
    short_list = [item for item in short_list if 'functional' not in item['field']]
    
    out_dir = Ldir['data'] / 'varinfo'
    Lfun.make_dir(out_dir)
    out_fn = out_dir / 'varinfo_list.p'
    pickle.dump(short_list, open(out_fn, 'wb'))
    
def get_varinfo(vn, vartype='state'):
    """
    This looks through the pre-parsed varinfo.yaml and returns a dict for a given variable.
    """
    Ldir = Lfun.Lstart()
    in_fn = Ldir['data'] / 'varinfo' / 'varinfo_list.p'
    if in_fn.is_file():
        short_list = pickle.load(open(in_fn, 'rb'))
    else:
        make_varinfo_list()
        short_list = pickle.load(open(in_fn, 'rb'))
        
    # get the dict for the requested variable
    short_list = [item for item in short_list if item['variable']==vn]
        
    # Associate grid_type with a tuple of spatial dimensions
    grid_type_dict = {
            'r2dvar': ('eta_rho', 'xi_rho'),
            'u2dvar': ('eta_u', 'xi_u'),
            'v2dvar': ('eta_v', 'xi_v'),
            'r3dvar': ('s_rho', 'eta_rho', 'xi_rho'),
            'u3dvar': ('s_rho', 'eta_u', 'xi_u'),
            'v3dvar': ('s_rho', 'eta_v', 'xi_v'),
            'nulvar': ('BLANK',)
            }
    
    if vartype=='state':
        short_list = [item for item in short_list if 'climatology' not in item['long_name']]
    elif vartype=='climatology':
        if vn in ['zeta', 'u', 'v', 'ubar', 'vbar']:
            # only these have explicity climatology versions
            short_list = [item for item in short_list if 'climatology' in item['long_name']]
        else:
            # the rest are the same as state variables
            short_list = [item for item in short_list if 'climatology' not in item['long_name']]
    else:
        print('Error in zrfun.get_varinfo(), unknown vartype: ' + vartype)
        sys.exit()
    
    if len(short_list)==1:
        vinfo = short_list[0]
        vinfo['space_dims_tup'] = grid_type_dict[vinfo['type']]
        vinfo['time_name'] = vinfo['time'] # for compatibility with forcing code
        
        # rename the time variable to work with climatology
        # NOTE: this does not work for the bry files - for these you just use vinfo['time'].
        if vartype=='climatology':
            if vn in ['zeta', 'u', 'v', 'ubar', 'vbar']:
                pass
            else:
                vinfo['time_name'] = vn + '_time'
        
        if vn in ['ubar', 'vbar']:
            vinfo['long_name'] = vinfo['long_name'].replace('integrated', 'averaged')
        return vinfo
    else:
        print('Error in zrfun.get_varinfo: vn not unique for ' + vn)
        vinfo = short_list
        return vinfo
    
"""
This is a dict to use for compression when saving an xarray Dataset, e.g. with lines like:
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    ds.to_netcdf(out_fn, encoding=Enc_dict)
Using compression (zlib=True, complevel=1) results in files that are just 2% of the
uncompressed files (for hc0, which has a lot of nan's).
Using complevel=9 makes the files half as big as complevel=1, but takes about 10x longer.
"""
enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}



