"""
Variables and functions used for plotting the multi-layer TEF extractions,
and for the flux code.

"""
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan

def get_two_layer(in_dir, sect_name):
    """
    Form time series of 2-layer TEF quantities, from the multi-layer bulk values.
    """
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    
    # determine which variables to process
    vn_list = []
    vec_list = []
    for vn in bulk.data_vars:
        if ('time' in bulk[vn].coords) and ('layer' in bulk[vn].coords):
            vn_list.append(vn)
        elif ('time' in bulk[vn].coords) and ('layer'  not in bulk[vn].coords):
            vec_list.append(vn)
    vn_list.remove('q') # transport is handled separately from tracers
    
    # separate positive and negative transports
    QQ = bulk.q.values
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan
    
    # form two-layer versions of volume and tracer transports
    Qp = np.nansum(QQp, axis=1)
    Qm = np.nansum(QQm, axis=1)
    if True:
        # Mask out times when the transport is too small to
        # use for tracer averaging.
        # Not needed? Seems like a possible source of budget errors.
        Qp[Qp<np.nanmean(Qp)/100] = np.nan
        Qm[Qm>np.nanmean(Qm)/100] = np.nan
    QCp = dict()
    QCm = dict()
    for vn in vn_list:
        QCp[vn] = np.nansum(QQp*(bulk[vn].values), axis=1)
        QCm[vn] = np.nansum(QQm*(bulk[vn].values), axis=1)
    
    # form flux-weighted tracer concentrations
    Cp = dict()
    Cm = dict()
    for vn in vn_list:
        Cp[vn] = QCp[vn]/Qp
        Cm[vn] = QCm[vn]/Qm
        
    # pack results in a DataFrame
    tef_df = pd.DataFrame(index=bulk.time.values)
    tef_df['q_p']=Qp
    tef_df['q_m']=Qm
    for vn in vn_list:
        tef_df[vn+'_p'] = Cp[vn]
        tef_df[vn+'_m'] = Cm[vn]
    # also pass back time series like qprism, for convenience
    for vn in vec_list:
        tef_df[vn] = bulk[vn].values
        
    bulk.close()
            
    return tef_df, vn_list, vec_list
    
def get_two_layer_from_list(bulk_dir, sntup_list):
    """
    This combines the transport of several sections in a consistent way.
    
    Assumes sntup_list is a list of tuples, e.g. [('ai1',-1),('dp',-1)].
    """
    for sntup in sntup_list:
        if not isinstance(sntup,tuple):
            print('Error from tef_fun.get_two_layer_from_list()')
            print('Need to pass a list of sns tulples.')
            sys.exit()
            
        sign_dict = dict()
        tef_df_dict = dict()
        
        vn_list = ['salt'] # NOTE: this will need to be generalized to more tracers!
        for sn_tup in sect_name:
            sn = sn_tup[0]
            sign_dict[sn] = sn_tup[1]
            tef_df_dict[sn], vn_list, vec_list = get_two_layer(bulk_dir, sn)
            
        ii = 1
        nsect = len(sect_name)
        for sn in tef_df_dict.keys():
            sgn = sign_dict[sn]
            if ii == 1:
                tef_df = tef_df_dict[sn].copy()
                tef_df[pd.isnull(tef_df)] = 0
                tef_df1 = tef_df.copy()
                for vn in vn_list:
                    if sgn == 1:
                        tef_df[vn+'_q_p'] = tef_df1['q_p'] * tef_df1[vn+'_p']
                        tef_df[vn+'_q_m'] = tef_df1['q_m'] * tef_df1[vn+'_m']
                    elif sgn == -1:
                        tef_df['q_p'] = -tef_df1['q_m']
                        tef_df['q_m'] = -tef_df1['q_p']
                        tef_df[vn+'_q_p'] = -tef_df1['q_m'] * tef_df1[vn+'_m']
                        tef_df[vn+'_q_m'] = -tef_df1['q_p'] * tef_df1[vn+'_p']
                tef_df['ssh'] *= 1/nsect
            else:
                tef_df1 = tef_df_dict[sn].copy()
                tef_df1[pd.isnull(tef_df1)] = 0
                for vn in vn_list:
                    if sgn == 1:
                        tef_df['q_p'] += tef_df1['q_p']
                        tef_df['q_m'] += tef_df1['q_m']
                        tef_df[vn+'_q_p'] += tef_df1['q_p'] * tef_df1[vn+'_p']
                        tef_df[vn+'_q_m'] += tef_df1['q_m'] * tef_df1[vn+'_m']
                    elif sgn == -1:
                        tef_df['q_p'] += -tef_df1['q_m']
                        tef_df['q_m'] += -tef_df1['q_p']
                        tef_df[vn+'_q_p'] += -tef_df1['q_m'] * tef_df1[vn+'_m']
                        tef_df[vn+'_q_m'] += -tef_df1['q_p'] * tef_df1[vn+'_p']
                for vn in ['qprism', 'qnet', 'fnet']:
                    tef_df[vn] += sgn * tef_df1[vn]
                tef_df['ssh'] += tef_df1['ssh']/nsect
            ii+= 1
        for vn in vn_list:
            tef_df[vn+'_p'] = tef_df[vn+'_q_p'] / tef_df['q_p']
            tef_df[vn+'_m'] = tef_df[vn+'_q_m'] / tef_df['q_m']
            
    return tef_df

# colors to associate with each channel (the keys in channel_ and seg_dict)
clist = ['blue', 'red', 'olive', 'orange']

units_dict = {'salt':'g/kg',
        'temp':'degC',
        'oxygen':'uM DO',
        'NO3':'uM N',
        'phytoplankton':'uM N',
        'zooplankton':'uM N',
        'detritus':'uM N',
        'Ldetritus':'uM N',
        'Ntot':'uM N',
        'TIC':'uM C',
        'alkalinity':'uM Eq'}

def make_dist(x,y):
    NS = len(x)
    xs = np.zeros(NS)
    ys = np.zeros(NS)
    xs, ys = zfun.ll2xy(x, y, x[0], y[0])
    dx = np.diff(xs)
    dy = np.diff(ys)
    dd = (dx**2 + dy**2)**.5 # not clear why np.sqrt throws an error
    dist = np.zeros(NS)
    dist[1:] = np.cumsum(dd/1000) # convert m to km
    return dist
        