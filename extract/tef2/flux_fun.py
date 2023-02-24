"""
Variables and functions used for plotting the multi-layer TEF extractions,
and for the flux code.

"""
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import pickle

from lo_tools import Lfun, zfun
# import tef_fun

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan

def get_two_layer(in_dir, sect_name):
    """
    Form time series of 2-layer TEF quantities, from the multi-layer bulk values.
    """
    
    bulk = pickle.load(open(in_dir / sect_name, 'rb'))
        
        
    QQ = bulk['q']
    dt = bulk['ot'] # datetimes
    dti = pd.DatetimeIndex(dt)
    
    vn_list = ['salt'] # need to fix this so it handles all variables automatically
        
    # separate positive and negative transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan
    
    # form two-layer versions of volume and tracer transports
    Qp = np.nansum(QQp, axis=1)
    Qm = np.nansum(QQm, axis=1)
    # mask out times when the transport is too small to use for tracer averaging
    Qp[Qp<np.nanmean(Qp)/10] = np.nan
    Qm[Qm>np.nanmean(Qm)/10] = np.nan
    QCp = dict()
    QCm = dict()
    for vn in vn_list:
        QCp[vn] = np.nansum(QQp*(bulk[vn]), axis=1)
        QCm[vn] = np.nansum(QQm*(bulk[vn]), axis=1)
    
    # form flux-weighted tracer concentrations
    Cp = dict()
    Cm = dict()
    for vn in vn_list:
        Cp[vn] = QCp[vn]/Qp
        Cm[vn] = QCm[vn]/Qm
    
    tef_df = pd.DataFrame(index=dti)
    tef_df['q_p']=Qp
    tef_df['q_m']=Qm
    tef_df['qabs'] = bulk['qabs']
    tef_df['qnet'] = bulk['qnet']
    tef_df['fnet'] = bulk['fnet']
    tef_df['ssh'] = bulk['ssh']
    
    for vn in vn_list:
        tef_df[vn+'_p'] = Cp[vn]
        tef_df[vn+'_m'] = Cm[vn]
            
    return tef_df

            
# colors to associate with each channel (the keys in channel_ and seg_dict)
clist = ['blue', 'red', 'olive', 'orange']

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
        