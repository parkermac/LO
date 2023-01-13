"""
Module of functions for LO/obs and associated code.
"""
import numpy as np
import pandas as pd

def renumber_cid(df):
    # Rework cid (cast ID) to be increasing from zero in steps of one.
    a = df.cid.values
    au = df.cid.unique() # returns uniques in order
    u_dict = dict(zip(au, np.arange(len(au))))
    b = np.nan * np.ones(len(a))
    for ii in u_dict.keys():
        b[a==ii] = u_dict[ii]
    df['cid'] = b
    
    return df
    
def make_info_df(df):
    # Also pull out a dateframe with station info to use for model cast extractions.
    ind = df.cid.unique()
    col_list = ['lon','lat','time','name','cruise']
    info_df = pd.DataFrame(index=ind, columns=col_list)
    for cid in df.cid.unique():
        info_df.loc[cid,col_list] = df.loc[df.cid==cid,col_list].iloc[0,:]
    info_df.index.name = 'cid'
    info_df['time'] = pd.to_datetime(info_df['time'])
    
    return info_df