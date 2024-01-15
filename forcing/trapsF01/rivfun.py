"""
Some extra river functions.

"""

import sys
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from lo_tools import Lfun
from lo_tools import river_functions as rivf

def get_tc_rn(ri_df):
    """
    Makes a new column in df called 'tc_rn' which is the name of a
    river that has T climatology.
    
    As of 2021.04.19 we have a number of new rivers with T climatology,
    so this should be updated.  However it will work as is.
    """
    for rn in ri_df.index:
        if rn in ['coquille']:
            tc_rn = 'umpqua'
        elif rn in ['wilson']:
            tc_rn = 'nehalem'
        elif rn in ['naselle', 'willapa', 'chehalis', 'humptulips',
                    'queets', 'hoh', 'calawah', 'hoko', 'elwha',
                    'dungeness']:
            tc_rn = 'quinault'
        elif rn in ['dosewallips', 'duckabush', 'hamma', 'skokomish', 'deschutes',
                    'nisqually', 'puyallup', 'green', 'snohomish',
                    'stillaguamish']:#, 'nf_skokomish', 'sf_skokomish']:
            tc_rn = 'cedar'
        elif rn in ['samish']:
            tc_rn = 'nooksack'
        elif rn in ['squamish']:
            tc_rn = 'clowhom'
        elif rn in ['oyster', 'tsolum', 'englishman', 'cowichan']:
            tc_rn = 'nanaimo'
        elif rn in ['gold']:
            tc_rn = 'sarita'
        else:
            tc_rn = rn
        ri_df.loc[rn, 'tc_rn'] = tc_rn
    return ri_df
    
def get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days):
    # load historical and climatological data
    Hflow_df = pd.read_pickle(Ldir['Hflow_fn'])
    Cflow_df = pd.read_pickle(Ldir['Cflow_fn'])
    Ctemp_df = pd.read_pickle(Ldir['Ctemp_fn'])
    
    # initialize output dict
    qt_df_dict = dict()
    for rn in gri_df.index:
        
        # initialize screen output
        sout = '-- ' + rn + ': '
        
        rs = ri_df.loc[rn].copy() # a series with info for this river
        rs['got_data'] = False
        # initialize a qt (flow and temperature vs. time) DataFrame for this river
        qt_df = pd.DataFrame(index=dt_ind, columns=['clim','his','usgs','ec','nws','final','temperature'])
        # fill with historical and climatological fields
        qt_df.loc[:, 'clim'] = Cflow_df.loc[yd_ind,rn].values
        qt_df.loc[:, 'temperature'] = Ctemp_df.loc[yd_ind,ri_df.loc[rn,'tc_rn']].values
        try:
            qt_df.loc[:, 'his'] = Hflow_df.loc[dt_ind,rn].copy()
        except KeyError:
            pass # recent times will not be in the historical record

        if pd.notnull(qt_df.loc[:, 'his']).all():
            qt_df['final'] = qt_df['his']
            sout += 'filled from historical '
            rs.got_data = True
            
        # otherwise try (sequentially) to fill from
        # nws, or usgs, or ec
        if pd.notnull(rs.nws) and (Ldir['run_type'] == 'forecast'):
            rs, qt = rivf.get_nws_data(rs)
            if rs.got_data:
                qt_df['nws'] = qt.reindex(dt_ind)
                qt_df['final'] = qt_df['nws']
                sout += 'filled from nws forecast '
                
        if (not rs.got_data) and pd.notnull(rs.usgs):
            if rn in ['skokomish', 'hamma']:
                rs, qt = rivf.get_usgs_data_custom(rs, days)
            else:
                rs, qt = rivf.get_usgs_data(rs, days)
            if rs.got_data:
                qt_df['usgs'] = qt
                qt_df['final'] = qt_df['usgs']
                sout += 'filled from usgs '
                
        if (not rs.got_data) and pd.notnull(rs.ec):
            rs, qt = rivf.get_ec_data(rs, days)
            if rs.got_data:
                qt_df['ec'] = qt
                qt_df['final'] = qt_df['ec']
                sout += 'filled from ec '
        
        # check results and fill with extrapolation (ffill) or climatology
        if ( pd.isnull(qt_df['final'].values).any() and
                not pd.isnull(qt_df['final'].values).all() ):
            qt_df['final'] = qt_df['final'].ffill(axis=0)
            sout += ':: extended by ffill'
        if pd.isnull(qt_df['final'].values).any():
            qt_df['final'] = qt_df['clim']
            sout += '!! still missing values after ffill - all filled with climatology'
        if (qt_df['final'].values < 0).any():
            qt_df['final'] = qt_df['clim']
            sout += '!! negative values - all filled with climatology'
        if pd.isnull(qt_df['final'].values).any():
            sout += '>>>>>>> flow has missing values!! <<<<<<<<<'
        # Temperature data
        if pd.isnull(qt_df['temperature'].values).any():
            sout += '>>>>>>> temp has missing values!! <<<<<<<<<'
        # save in the dict
        qt_df_dict[rn] = qt_df
        
        # screen output
        print(sout)
        sys.stdout.flush()
        
    return qt_df_dict
    
def get_bio_vec(vn, rn, yd_ind):
    ndt = len(yd_ind)
    yd = yd_ind.values
    ovec = np.ones(ndt)
    if vn == 'NO3':
        if rn == 'fraser':
            vv = 2 + (13/2) + (13/2)*np.cos(2*np.pi*((yd-30)/366))
        elif rn == 'columbia':
            vv = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
        else:
            vv = 5 * ovec
    elif vn == 'Oxyg':
        vv = 350 * ovec
    elif vn in ['TAlk', 'TIC']:
        if rn in ['columbia', 'deschutes', 'duwamish']:
            vv = 1000 * ovec
        else:
            vv = 300 * ovec
    else:
        vv = 0 * ovec # all others filled with zeros
    return vv    
    
