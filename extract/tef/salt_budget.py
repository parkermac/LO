"""
Do volume and salt budgets for user-specified volumes.

Run with a command like:
run salt_budget -test True

"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))

import Lfun
import tef_fun
import flux_fun

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta
import argparse
import xarray as xr

# debugging imports
import netCDF4 as nc
from time import time
import zfun

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-test', '--testing', type=zfun.boolean_string, default=False)
args = parser.parse_args()
testing = args.testing

# Get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)

if testing:
    vol_list = ['Puget Sound']
    save_figs = False
else:
    vol_list = ['Salish Sea', 'Puget Sound', 'Hood Canal']
    save_figs = True

plt.close('all')
year = 2018
for which_vol in vol_list:

    year_str = str(year)
    date_str = '_' + year_str + '.01.01_' + year_str + '.12.31'
    
    riv_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_roms' / ('extraction' + date_str + '.p')
    tef_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / ('bulk' + date_str)
    seg_fn = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / ('segments_ds' + date_str + '.nc')
    vol_dir = Ldir['LOo'] / 'extract' / 'tef' / ('volumes_' + Ldir['gridname'])

    # Info specific to each volume
    # The sign for each section indicates which direction is INTO the volume.
    if which_vol == 'Salish Sea':
        seg_list = list(v_lp_df.columns)
        sect_sign_dict = {'jdf1':1, 'sog5':-1}
    elif which_vol == 'Puget Sound':
        seg_list = (flux_fun.ssA + flux_fun.ssM + flux_fun.ssT
            + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH)
        sect_sign_dict = {'ai1':1, 'dp':1}
    elif which_vol == 'Hood Canal':
        seg_list = flux_fun.ssH
        sect_sign_dict = {'hc1':-1}

    # SECTION INFO
    sect_df = tef_fun.get_sect_df(Ldir['gridname'])

    # RIVERS
    river_list = []
    for seg_name in seg_list:
        seg = flux_fun.segs[seg_name]
        river_list = river_list + seg['R']
    riv_df = pd.read_pickle(riv_fn)
    riv_df.index += timedelta(days=0.5)
    riv_df = riv_df[river_list]

    # TEF at SECTIONS
    tef_df_dict = dict()
    sect_list = list(sect_sign_dict.keys())
    for sn in sect_list:
        tef_df_dict[sn], in_sign, _, _ = flux_fun.get_two_layer(tef_dir, sn, Ldir['gridname'])
        if in_sign != sect_sign_dict[sn]:
            print('WARNING: potential sign error!!')
            
    # SEGMENT TIME SERIES
    pad = 36
    
    seg_ds = xr.open_dataset(seg_fn)
    seg_ds = seg_ds.sel(seg=seg_list)
    
    seg_NT = len(seg_ds.coords['time'])
    nanvec = np.nan * np.ones(seg_NT)
    
    vt = nanvec.copy()
    vt[1:-1] = (seg_ds.volume[2:].values - seg_ds.volume[:-2].values).sum(axis=1)/(2*3600)
    vt_lp = zfun.lowpass(vt, f='godin')[pad:-pad+1:24]
    
    seg_vn_list = list(seg_ds.data_vars.keys())
    seg_vn_list.remove('volume')
    
    cvt_lp_dict = {}
    for vn in seg_vn_list:
        cvt = nanvec.copy()
        cvt[1:-1] = (seg_ds.volume[2:].values*seg_ds[vn][2:].values
            - seg_ds.volume[:-2].values*seg_ds[vn][:-2].values).sum(axis=1)/(2*3600)
        cvt_lp_dict[vn] = zfun.lowpass(cvt, f='godin')[pad:-pad+1:24]
                        
    # BUDGETS
    indall = tef_df_dict[sect_list[0]].index

    # volume budget
    vol_df = pd.DataFrame(0, index=indall, columns=['Qin','Qout'])
    for sect_name in sect_list:
        df = tef_df_dict[sect_name]
        vol_df['Qin'] = vol_df['Qin'] + df['Qin']
        vol_df['Qout'] = vol_df['Qout'] + df['Qout']
    vol_df['Qr'] = riv_df.sum(axis=1)
    vol_df.loc[:, 'dV_dt'] = vt_lp
    vol_df['Error'] = vol_df['dV_dt'] - vol_df.loc[:,'Qin'] - vol_df.loc[:,'Qout'] - vol_df.loc[:,'Qr']
    vol_rel_err = vol_df['Error'].mean()/vol_df['Qr'].mean()

    # tracer budgets
    C = dict()
    for vn in seg_vn_list:
        c_df = pd.DataFrame(0, index=indall, columns=['QCnet'])#['QCin','QCout'])
        for sect_name in sect_list:
            df = tef_df_dict[sect_name]
            c_df['QCnet'] = c_df['QCnet'] + df['Qin']*df[vn+'_in'] + df['Qout']*df[vn+'_out']
            # c_df['QCin'] = c_df['QCin'] + df['Qin']*df[vn+'_in']
            # c_df['QCout'] = c_df['QCout'] + df['Qout']*df[vn+'_out']
        c_df.loc[:, 'dCnet_dt'] = cvt_lp_dict[vn]
        c_df['Source'] = c_df['dCnet_dt'] - c_df['QCnet']
        #c_df['Source'] = c_df['dCnet_dt'] - c_df['QCin'] - c_df['QCout']
        C[vn] = c_df
        
    C['Ntot'] = C['NO3']+C['phytoplankton']+C['zooplankton']+C['detritus']+C['Ldetritus']

