"""
Module of functions for LO/obs and associated code.
"""
import numpy as np

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