"""
Code to convert bottle obsmod data to json format.

This relies on dicts of dataframes that have been previously
created by LO/obsmod/combin_obs_mod.py.
"""

from lo_tools import Lfun
import pandas as pd
import numpy as np
from lo_tools import obs_functions

Ldir = Lfun.Lstart()

# input and output locations
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'

out_dir = Ldir['LOo'] / 'obsmod_json'
Lfun.make_dir(out_dir)

# load data
year_list = range(2013,2026)
#year_list = [2014]

source = 'combined'
otype = 'bottle'
gtagex = 'cas7_t1_x11ab'

for year in year_list:

    print(year)

    # load data
    in_fn = in_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '.p')
    df_dict = pd.read_pickle(in_fn)

    # save parts to json

    # data columns to include
    info_list = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z', 'source']
    fld_list = ['CT', 'SA', 'DO (uM)', 'NO3 (uM)', 'DIC (uM)', 'TA (uM)']
    full_list = info_list + fld_list

    # region to include (cas7 model domain)
    aa = [-130, -122, 42, 52]
    x0, x1, y0, y1 = aa

    df_dict2 = {}
    for om in ['obs', 'mod']:
        if om == 'obs':
            gtx = 'obs'
        elif om == 'mod':
            gtx = gtagex
        # select Dataframe
        df = df_dict[gtx].copy()
        # keep only data in the LO domain
        df = df[(df['lon']>x0) & (df['lon']<x1) & (df['lat']>y0) & (df['lat']<y1)]
        # make sure we have all required columns, even if there is no data
        for fld in fld_list:
            if fld not in df.columns:
                df.loc[:,fld] = np.nan
        # keep only desired columns
        df = df[full_list]
        # and save for subsequent processing
        df_dict2[om] = df

    dfo = df_dict2['obs'].copy()
    dfm = df_dict2['mod'].copy()

    nfld = len(fld_list)
    obad = dfo[fld_list].isnull().sum(axis=1) == nfld
    mbad = dfm[fld_list].isnull().sum(axis=1) == nfld
    allgood = ~obad & ~mbad
    dfo = dfo.loc[allgood,:]
    dfm = dfm.loc[allgood,:]
    df_dict3 = {}
    df_dict3['obs'] = dfo
    df_dict3['mod'] = dfm
    print(' dfo length = ' + str(len(dfo)))

    for om in ['obs', 'mod']:
        out_fn = out_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '_' + om + '.json')
        df_dict3[om].to_json(out_fn, double_precision=3)

        # also pull out station locations with cid as index
        if om == 'obs':
            info_df = obs_functions.make_info_df(df)
            info_df['cid'] = info_df.index.to_numpy()
            info_df.index = info_df.index.astype(int)
            out_info_fn = out_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '_info.json')
            info_df.to_json(out_info_fn, double_precision=3)

    # obs_df = df_dict['obs']
    # # keep only data in the LO domain
    # obs_df = obs_df[(df['lon']>x0) & (df['lon']<x1) & (df['lat']>y0) & (df['lat']<y1)]
    # # make sure we have all required columns, even if there is no data
    # for fld in fld_list:
    #     if fld not in obs_df.columns:
    #         obs_df[fld] = np.NaN
    # obs_df = obs_df[full_list]
    # out_obs_fn = out_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '_obs.json')
    # obs_df.to_json(out_obs_fn, double_precision=3)

    # mod_df = df_dict[gtagex]
    # for fld in fld_list:
    #     if fld not in mod_df.columns:
    #         mod_df[fld] = np.NaN
    # mod_df = mod_df[full_list]
    # out_mod_fn = out_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '_mod.json')
    # mod_df.to_json(out_mod_fn, double_precision=3)

    # # also pull out station locations with cid as index
    # info_df = obs_functions.make_info_df(obs_df)
    # info_df['cid'] = info_df.index.to_numpy()
    # info_df.index = info_df.index.astype(int)
    # out_info_fn = out_dir / (source + '_' + otype + '_' + str(year) + '_' + gtagex + '_info.json')
    # info_df.to_json(out_info_fn, double_precision=3)

    # # output looks like:
    # # 
