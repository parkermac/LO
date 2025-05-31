"""
Helper funtions to generate traps forcing
"""

import numpy as np
import pandas as pd
import sys
from lo_tools import forcing_argfun2 as ffun


Ldir = ffun.intro() # this handles all the argument passing

def get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type, trapsD):

    # Only add biology to pre-existing LO river if Ecology has data
    if traps_type == 'LOriv':
        # get names of duplicate rivers
        repeatrivs_fn = Ldir['data'] / trapsD / 'LiveOcean_SSM_rivers.xlsx'
        repeatrivs_df = pd.read_excel(repeatrivs_fn)
        LObio_names_all = list(repeatrivs_df.loc[repeatrivs_df['in_both'] == 1, 'LO_rname'])
        # remove the weird rivers
        weird_duplicate_rivers = ['Alberni Inlet', 'Chehalis R', 'Gold River', 'Willapa R', 'Columbia R', 'Comox']
        # Note that these are the names that LO calls the rivers
        LObio_names = [rname for rname in LObio_names_all if LO2SSM_name(rname,trapsD) not in weird_duplicate_rivers]

    # load climatological data
    if traps_type != 'LOriv':
        # don't need flow and temp for pre-existing LO rivers
        Cflow_df = pd.read_pickle(Ldir['Cflow_'+traps_type+'_fn'])
        Ctemp_df = pd.read_pickle(Ldir['Ctemp_'+traps_type+'_fn'])
    CDO_df   = pd.read_pickle(Ldir['CDO_'+traps_type+'_fn'])
    CNH4_df  = pd.read_pickle(Ldir['CNH4_'+traps_type+'_fn'])
    CNO3_df  = pd.read_pickle(Ldir['CNO3_'+traps_type+'_fn'])
    CTalk_df = pd.read_pickle(Ldir['CTalk_'+traps_type+'_fn'])
    CTIC_df  = pd.read_pickle(Ldir['CTIC_'+traps_type+'_fn'])

    # year day index starts from 1. Convert to start from 0 to work with python
    yd_ind = yd_ind - 1

    # initialize output dict
    qtbio_df_dict = dict()
    for rn in gri_df.index:

        # convert LO river name to SSM river name
        if traps_type == 'LOriv':
            if rn in LObio_names:
                rn = LO2SSM_name(rn,trapsD)
            else:
                # skips rivers for which Ecology does not have data
                continue    
        
        # initialize screen output
        sout = '-- ' + str(rn) + ': '
        
        # initialize a qtbio (flow and temperature vs. time) DataFrame for this river
        qtbio_df = pd.DataFrame(index=dt_ind, columns=['flow','temp','Oxyg','NH4','NO3','TAlk','TIC'])
        # fill with climatological fields
        if traps_type != 'LOriv':
            qtbio_df.loc[:, 'flow'] = Cflow_df.loc[yd_ind,rn].values
            qtbio_df.loc[:, 'temp'] = Ctemp_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'Oxyg']   = CDO_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'NH4']  = CNH4_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'NO3']  = CNO3_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'TAlk'] = CTalk_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'TIC']  = CTIC_df.loc[yd_ind,rn].values
        sout += 'filled from climatology (Mohamedali et al., 2020 dataset)'
          
        # save in the dict
        qtbio_df_dict[rn] = qtbio_df
        
        # screen output
        print(sout)
        sys.stdout.flush()
        
    return qtbio_df_dict

def combine_adjacent(lst):
    """
    Given a list, e.g. ['a','b','c','d']
    returns: ['a+b', 'c+d']
    """
    combined = [x + '+' + y for x, y in zip(lst[::2],lst[1::2])]
    return combined

def weighted_average(vn,qtbio_df_1, qtbio_df_2):
    '''
    Calculate the weighted average properties based on flowrate of two overlapping sources
    '''
    # get flowrates
    flow1 = qtbio_df_1['flow'].values
    flow2 = qtbio_df_2['flow'].values
    # get variable
    var1 = qtbio_df_1[vn].values
    var2 = qtbio_df_2[vn].values
    # calculate weighted average based on flowrate
    waverage = [np.average([var1[i], var2[i]], weights = [flow1[i], flow2[i]]) for i in range(len(flow1))]
    return waverage

def LO2SSM_name(rname,trapsD):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / trapsD / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]

    return rname_SSM