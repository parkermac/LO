"""
Do volume and tracer budgets for user-specified volumes.

To test on mac:
run tracer_budget -gtx cas7_trapsV00_meV00 -ctag c0 -riv trapsV00 -0 2017.07.04 -1 2017.07.06 -test True
run tracer_budget -gtx cas7_trapsV00_meV00 -ctag c0 -riv trapsV00 -0 2017.01.01 -1 2017.01.10 -test True
run tracer_budget -gtx cas7_trapsV00_meV00 -ctag c0 -riv trapsV00 -0 2017.01.01 -1 2017.12.31 -test True

Run with -test True to make a plot, and avoid overwriting the output directory.

Performance: fast.

"""

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun

import numpy as np
import pickle
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta

# debugging imports
from time import time

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

testing = Ldir['testing']
sect_gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
riv_gctag = Ldir['gridname'] + '_' + Ldir['riv']
date_str = '_' + Ldir['ds0'] + '_' + Ldir['ds1']

# get the budget_functions module, looking first in LO_user
pth = Ldir['LO'] / 'extract' / 'tef2'
upth = Ldir['LOu'] / 'extract' / 'tef2'
if (upth / 'budget_functions.py').is_file():
    print('Importing budget_functions from LO_user')
    bfun = Lfun.module_from_file('budget_functions', upth / 'budget_functions.py')
else:
    print('Importing budget_functions from LO')
    bfun = Lfun.module_from_file('budget_functions', pth / 'budget_functions.py')

if testing:
    vol_list = ['Puget Sound']
else:
    vol_list = ['Salish Sea', 'Puget Sound', 'Hood Canal']
    
if testing:
    import matplotlib.pyplot as plt
    from lo_tools import plotting_functions as pfun
    plt.close('all')

for which_vol in vol_list:
    vol_str = which_vol.replace(' ','_')

    year_str = str(datetime.strptime(Ldir['ds0'],Lfun.ds_fmt).year)

    # base location for output
    dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
        
    # GATHERING INPUT FILES
    
    # we use the "bulk" output to get transport through the open boundaries
    bulk_dir = dir0 / ('bulk' + date_str)
    
    # we use the "segments" output to get the net tracer storage in the volume
    seg_ds_fn = dir0 / ('segments' + date_str + '_' + sect_gctag + '_' + Ldir['riv'] + '.nc')
    seg_ds = xr.open_dataset(seg_ds_fn)
    
    # we use the river extraction to get river contributions to the budget
    dir1 = Ldir['LOo'] / 'pre' / 'river1'
    riv_ds_fn = dir1 / riv_gctag / 'Data_roms' / ('extraction' + date_str + '.nc')
    riv_ds = xr.open_dataset(riv_ds_fn)
    
     # we need the seg_info_dict to figure out which rivers to include, and which segments
    dir2 = Ldir['LOo'] / 'extract' / 'tef2'
    seg_info_dict_fn = dir2 / ('seg_info_dict_' + sect_gctag + '_' + Ldir['riv'] + '.p')
    seg_info_dict = pd.read_pickle(seg_info_dict_fn)
    
    # get the sect df to find all the valid section names
    sect_df_fn = dir2 / ('sect_df_' + sect_gctag + '.p')
    sect_df = pd.read_pickle(sect_df_fn)
    sn_list = list(sect_df.sn)
    
    # get info about bounding  and interior sections
    sntup_list, sect_base_list, outer_sns_list = bfun.get_sntup_list(sect_gctag, which_vol)
    
    # FIND WHICH SEGMENTS ARE PART OF THE VOLUME, AND WHICH RIVERS
    
    # start by figuring out all valid sns
    sns_list = []
    for snb in sect_base_list:
        for sn in sn_list:
            if snb in sn:
                for pm in ['_p','_m']:
                    sns = sn + pm
                    if (sns not in outer_sns_list) and (sns not in sns_list):
                        sns_list.append(sns)
                    else:
                        pass
                        # print(' - excluding ' + sns + ' from sns_list')
    # then use the valid sns to include only valid segments
    good_seg_key_list = []
    for sk in seg_info_dict.keys():
        this_sns_list = seg_info_dict[sk]['sns_list']
        check_list = [item for item in this_sns_list if item in sns_list]
        if len(check_list) >= 1:
            good_seg_key_list.append(sk)
        elif len(check_list) == 0:
            pass
    # and now that we have the segment keys we can also get the list of rivers
    good_riv_list = []
    for sk in good_seg_key_list:
        good_riv_list += seg_info_dict[sk]['riv_list']
        
    # ORGANIZE BUDGET TERMS
    
    vn_list = tef_fun.vn_list # get the default list
    seg_vn_list = [] # this will hold the ones we can make budgets for
    for vn in seg_ds.data_vars:
        if ('time' in seg_ds[vn].coords) and ('seg' in seg_ds[vn].coords):
            seg_vn_list.append(vn)
    # trim vn_list
    for vn in vn_list:
        if vn not in seg_vn_list:
            vn_list.remove(vn)
    if testing:
        # further trimming of vn_list
        vn_list = ['salt','temp']
    # add volume to vn_list
    vn_list += ['volume']
     
    # Set up a place to store things, and some helpful arrays:
    B_dict = dict() # a dict to hold separate budgets, each a DataFrame
    pad = 36 # for removing ends after Godin filter
    zvec_h = np.zeros(len(seg_ds.time)) # hourly
    nanvec_h = np.nan * zvec_h
    zvec_d = np.zeros(len(riv_ds.time)) # daily
    nanvec_d = np.nan * zvec_d
    # Note: we assume all the extractions cover the same time span, and that
    # seg_ds has the full list of hours, while riv_ds has the full list of days
    
    # CALCULATE BUDGETS FOR EACH VARIABLE
    
    for vn in vn_list:
        df = pd.DataFrame(0, index=riv_ds.time, columns=['net','d_dt','riv','ocn','surf','err'])
        
        # Rivers
        # These start as daily at noon, including the first and last days
        riv = zvec_d.copy()
        for rn in good_riv_list:
            if vn == 'volume':
                riv += riv_ds.sel(riv=rn).transport.values
            else:
                if vn in tef_fun.ocean_to_river_dict.keys():
                    vnr = tef_fun.ocean_to_river_dict[vn]
                else:
                    vnr = vn
                riv += riv_ds.sel(riv=rn).transport.values * riv_ds.sel(riv=rn)[vnr].values
        df['riv'] = riv
    
        # Segments
        # These start as hourly
        net_h = zvec_h.copy()
        surf_h = zvec_h.copy()
        for sk in good_seg_key_list:
            this_ds = seg_ds.sel(seg=sk)
            if vn == 'volume':
                net_h += this_ds.volume.values
            else:
                net_h += this_ds[vn].values * this_ds.volume.values
                if vn == 'salt':
                    # include EminusP
                    surf_h += this_ds.salt_surf.values * this_ds.area.values * this_ds.EminusP.values
                elif vn == 'temp':
                    # include shflux
                    rho = 1025 # approximate density [kg m-3]
                    Cp = 4000 # approximate heat capacity [J kg-1 degK-1]
                    surf_h += this_ds.area.values * this_ds.shflux.values / (rho * Cp)
                    # Units of surf_sh are [degC m3 s-1]
                    # Note: this does not work very well. May need to do a proper
                    # heat budget keeping track of actual Joules m-3, but this implies
                    # intervening in the extractions.
        d_dt_h= (net_h[2:] - net_h[:-2])/(2*3600)
        d_dt_h_full = nanvec_h.copy()
        d_dt_h_full[1:-1] = d_dt_h
        d_dt = nanvec_d.copy()
        d_dt[1:-1] = zfun.lowpass(d_dt_h_full, f='godin')[pad:-pad+1:24]
        df['d_dt'] = d_dt
        net = nanvec_d.copy()
        net[1:-1] = zfun.lowpass(net_h, f='godin')[pad:-pad+1:24]
        df['net'] = net
        surf = nanvec_d.copy()
        surf[1:-1] = zfun.lowpass(surf_h, f='godin')[pad:-pad+1:24]
        df['surf'] = surf
        
        # Ocean: transport through bounding sections
        # These start as low-passed, daily at noon, missing the first and last days
        ocn = zvec_d.copy()[1:-1]
        for tup in sntup_list: # add up contribution of all sections
            sn = tup[0]
            sgn = tup[1]
            bulk= xr.open_dataset(bulk_dir / (sn + '.nc'))
            qnet = bulk.qnet.values
            if vn == 'volume':
                ocn += qnet * sgn
            else:
                ocn += sgn * np.nansum((bulk.q * bulk[vn]).values, axis=1)
            bulk.close()
        ocn_full = nanvec_d.copy()
        ocn_full[1:-1] = ocn
        df['ocn'] = ocn_full
        
        # Error, or non-conservation
        df['err'] = df.d_dt - df.riv - df.ocn - df.surf
        
        # save the results for this variable to the output dict
        B_dict[vn] = df.copy()
        
    if testing:
        
        pfun.start_plot(figsize=(12,8))
        
        for vn in vn_list:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            this_df = B_dict[vn]
            this_df.loc[:,['d_dt','riv','ocn','surf','err']].plot(ax=ax, grid=True)
            ax.set_xlim(this_df.index[0],this_df.index[-1])
            ax.set_title(vn + ' budget: ' + which_vol)
            ax.set_ylabel(tef_fun.units_dict[vn] + ' m3 s-1')
            
        plt.show()
        pfun.end_plot()
        
    else:
        # save output as a pickled dict of DataFrames
        out_dir = dir0 / ('Budgets_' + date_str)
        Lfun.make_dir(out_dir, clean=True)
        pd.to_pickle(B_dict, (out_dir / (vol_str + '.p')))
        

