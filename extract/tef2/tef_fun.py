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

# Default list of tracers to get, for the new ROMS runs, as of 2023.
vn_list = ['salt', 'temp', 'oxygen','chlorophyll',
            'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN',
            'TIC', 'alkalinity', 'SdetritusC', 'LdetritusC']
            
# Dict that associates tracer names in the river forcing with those in
# the history files. Note that in the river forcing file they are all
# called river_[variable name], e.g. river_TAlk.
# This ONLY includes names that are DIFFERENT between river and history files.
river_to_ocean_dict = {
    'Phyt':'phytoplankton',
    'Zoop':'zooplankton',
    'LDeN':'LdetritusN',
    'SDeN':'SdetritusN',
    'Chlo':'chlorophyll',
    'TAlk':'alkalinity',
    'LDeC':'LdetritusC',
    'SDeC':'SdetritusC',
    'Oxyg':'oxygen'}
# same dict reversing keys and values
ocean_to_river_dict = {river_to_ocean_dict[k]:k for k in river_to_ocean_dict.keys()}

units_dict = {
    'volume':'',
    'salt':'g/kg',
    'temp':'degC',
    'oxygen':'uM DO',
    'NO3':'uM N',
    'phytoplankton':'uM N',
    'chlorophyll':'mg m-3',
    'zooplankton':'uM N',
    'SdetritusN':'uM N',
    'LdetritusN':'uM N',
    'SdetritusC':'uM C',
    'LdetritusC':'uM C',
    'Ntot':'uM N',
    'TIC':'uM C',
    'alkalinity':'uM Eq'}


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
        