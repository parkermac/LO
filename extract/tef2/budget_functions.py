"""
This is a module of functions used for budget calculations. It is a likely
place where users will need to make their own version, so in code that
uses it we add a "hook" to look for a version in LO_user.

"""
import flux_fun # how to make this work in LO_user?
import pandas as pd

def get_sect_list(gctag, vol_name):
    """
    This has lists of sections that are the open boundaries to selected
    volumes.
    
    NOTE: We handle combining sections by putting them in a tuple, and
    then each section is in its own tuple with a 1 or -1 to indicate the sign.
    The reason for the sign is that otherwise we have no way of knowing
    how to combine the two.
    """
    if gctag == 'cas6_c0':
        if vol_name == 'Puget Sound':
            sect_list = [(('ai1',-1),('dp',-1))] # the sign is for inflow direction
            bounding_sect_list = ['ai1_p','dp_p'] # exclude segments that have these
            seg_base_list = ['ai','wb','hc','mb','tn','ss']
    return sect_list, bounding_sect_list, seg_base_list

def get_two_layer_from_list(bulk_dir, sect_list):
    """
    The columns returned in the DataFrame from flux_fun.get_two_layer are:
    ['q_p', 'q_m', 'qprism', 'qnet', 'fnet', 'ssh', 'salt_p', 'salt_m']
    among others...
    """
    for sect_name in sect_list:
        if isinstance(sect_name,tuple):
            sign_dict = dict()
            tef_df_dict = dict()
            vn_list = ['salt'] # NOTE: this will need to be generalized to more tracers!
            for sn_tup in sect_name:
                sn = sn_tup[0]
                sign_dict[sn] = sn_tup[1]
                tef_df_dict[sn] = flux_fun.get_two_layer(bulk_dir, sn)
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
        else:
            tef_df = flux_fun.get_two_layer(in_dir, sect_name)
            
    return tef_df
    
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
