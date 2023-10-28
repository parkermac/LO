"""
Module of functions for LO/obs and associated code.
"""
import numpy as np
import pandas as pd

# Names of the seasons
season_name_dict = {0:'Winter',1:'Spring',2:'Summer',3:'Fall'}

# Colors for the seasons
season_c_dict = {0:'b',1:'g',2:'r',3:'orange'}

# dayofyear for the seasons (edges)
season_daylist = [0, 92,183,275,367]


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