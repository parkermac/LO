"""
Do volume and tracer budgets for user-specified volumes.

Run with a command like:
run tracer_budget -gtx cas6_v00_uu0m -ctag c0 -riv riv00 -0 2022.01.01 -1 2022.12.31 -test True
run tracer_budget -gtx cas7_trapsV00_meV00 -ctag c0 -riv trapsV00 -0 2017.01.01 -1 2017.01.10 -test True

"""

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun

import matplotlib.pyplot as plt
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

for which_vol in vol_list:
    vol_str = which_vol.replace(' ','_')

    year_str = str(datetime.strptime(Ldir['ds0'],Lfun.ds_fmt).year)

    # location for output
    dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
    if testing == False:
        out_dir = dir0 / ('Budgets_' + vol_str + date_str)
        Lfun.make_dir(out_dir, clean=True)
        
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
        
    # ORGANIZE BUDGET TERMS INTO A SINGLE PANDAS DATAFRAME
    
    # Set up some places to store things:
    # - from segments
    VC_dict = dict() # volume integrals vs. time, including d/dt of each, hourly
    VC_lp_dict = dict() # low-passed version of the items in VC_dict, daily
    # - from rivers, open boundary sections, and surface/bottom fluxes
    F_dict = dict()
    F_lp_dict = dict()
    
    # Use seg_ds to figure out which tracers we can process.
    vn_list = []
    for vn in seg_ds.data_vars:
        if ('time' in seg_ds[vn].coords) and ('seg' in seg_ds[vn].coords):
            vn_list.append(vn)
    vn_list.remove('volume') # we treat this separately from the other tracers
    
    # Rivers
    qr = riv_ds.sel(riv=good_riv_list).transport.sum('riv')
    # an xr DataArray with a time series (noon daily) of net river flow
    qr_ser = qr.to_series() # simple method for xr -> pd
    
    # Segment volume and dv/dt, and other tracers
    #
    vol_hourly = seg_ds.sel(seg=good_seg_key_list).volume.sum('seg') # a DataArray
    vol_vec = vol_hourly.values
    # rate of change of volume
    dvdt_vec = np.nan * vol_vec
    dvdt_vec[1:-1] = (vol_vec[2:] - vol_vec[:-2])/(2*3600)
    pad = 36
    vol_daily = zfun.lowpass(vol_vec, f='godin')[pad:-pad+1:24]
    dvdt_daily = zfun.lowpass(dvdt_vec, f='godin')[pad:-pad+1:24]
    vol_dt_daily = vol_hourly.time.values[pad:-pad+1:24]
    vol_ser = pd.Series(index=vol_dt_daily, data=vol_daily)
    dvdt_ser = pd.Series(index=vol_dt_daily, data=dvdt_daily)
    # practice keeping things organized
    VC_dict['volume'] = vol_vec
    VC_dict['dvdt'] = dvdt_vec
    VC_dict['time'] = seg_ds.time.values
    VC_lp_dict['volume'] = vol_daily
    VC_lp_dict['dvdt'] = dvdt_daily
    VC_lp_dict['time'] = vol_dt_daily
    # Also get the net tracer in all the segments. Note that what extract_segments.py
    # saves is the hourly average tracer in each segment.
    for vn in vn_list:
        for sk in good_seg_key_list:
            pass
    
    # Transport (low-passed, daily at noon)
    bulk_dir = dir0 / ('bulk' + date_str)
    ii = 0
    for tup in sntup_list: # add up contribution of all sections
        sn = tup[0]
        sgn = tup[1]
        bulk= xr.open_dataset(bulk_dir / (sn + '.nc'))
        qnet = bulk.qnet.values
        if ii == 0:
            qnet_vec = qnet * sgn
            otq = bulk.time.values
        else:
            qnet_vec += qnet * sgn
        ii += 1
    qnet_ser = pd.Series(index=otq, data=qnet_vec)
        
    # combine in a pandas DataFrame
    vol_df = pd.DataFrame()
    vol_df['riv'] = qr_ser # this has the longest time axis
    vol_df['vol'] = vol_ser
    vol_df['dvdt'] = dvdt_ser
    vol_df['qnet'] = qnet_ser
    vol_df['err'] = vol_df.dvdt - vol_df.riv - vol_df.qnet
    
    if testing:
        plt.close('all')
        vol_df.loc[:,['riv','dvdt','qnet','err']].plot()
        plt.show()


    #
    # # Tracer budgets
    # # F is the "flux" of a tracer, with units [tracer units]*m3/s
    # # Ftot, Fin, and Fout are at the ocen boundaries of the volume.  Fout is negative.
    # C = dict()
    # # The "normalized" budgets are averaged over a year, multiplied by a year of seconds,
    # # and divided by the mean volume, so they have units [tracer units].  So Cnorm['NO3']['dFnet_dt']
    # # is the change in total mean Nitrate over the year, and the other terms in Cnorm['NO3'] tell you where
    # # that change came from.  I'm not sure this is the right normalization to use.
    # Cnorm = dict()
    # for vn in tef_fun.vn_list:
    #     c_df = pd.DataFrame(0, index=indall, columns=['Ftot','Fin','Fout'])
    #     for sect_name in sect_list:
    #         df = tef_df_dict[sect_name]
    #         c_df['Ftot'] = c_df['Ftot'] + df['Qin']*df[vn+'_in'] + df['Qout']*df[vn+'_out']
    #         c_df['Fin'] = c_df['Fin'] + df['Qin']*df[vn+'_in']
    #         c_df['Fout'] = c_df['Fout'] + df['Qout']*df[vn+'_out']
    #     c_df['Fr'] = (riv_ds.transport * riv_ds[vn]).sum(axis=1)[1:-1]
    #     c_df.loc[:, 'dFnet_dt'] = cvt_lp_dict[vn]
    #     # the residual of the budget is assumed to be an unresolved Source or Sink (Sink is negative)
    #     # e.g. due to air-sea gas transfer, denitrification, or internal conversion to another tracer.
    #     c_df['Source/Sink'] = c_df['dFnet_dt'] - c_df['Ftot'] - c_df['Fr']
    #     C[vn] = c_df.copy()
    #     cn = c_df.mean()*(365*86400)/V
    #     cn = cn.rename({'dFnet_dt':'Change in Concentration', 'Ftot':'Inflow+Outflow',
    #         'Fin':'Inflow', 'Fout':'Outflow', 'Fr':'River'})
    #     cn = cn[['Change in Concentration', 'Inflow', 'Outflow', 'Inflow+Outflow', 'River','Source/Sink']]
    #     cn['Mean Concentration'] = vmean_dict[vn]
    #     Cnorm[vn] = cn
    #
    # C['Ntot'] = C['NO3']+C['phytoplankton']+C['zooplankton']+C['detritus']+C['Ldetritus']
    # Cnorm['Ntot'] = Cnorm['NO3']+Cnorm['phytoplankton']+Cnorm['zooplankton']+Cnorm['detritus']+Cnorm['Ldetritus']
    #
    # plt.close('all')
    # pfun.start_plot()
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # tstr = which_vol + ' Volume Budget [m3/s]'
    # vol_df.plot(ax=ax, grid=True, title=tstr)
    # if testing:
    #     plt.show()
    # else:
    #     fig.savefig(out_dir / 'volume.png')
    #
    # for vn in C.keys():
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     tstr = which_vol + ' ' + vn + ' Budget [' + units_dict[vn] + ' m3/s]'
    #     C[vn][['dFnet_dt','Ftot','Fr','Source/Sink']].plot(ax=ax, grid=True, title=tstr)
    #     if testing:
    #         pass
    #         #plt.show()
    #     else:
    #         fig.savefig(out_dir / (vn + '.png'))
    #
    # # text output
    # with open(out_dir / ('Annual_Mean_' + which_vol.replace(' ','_') + '_' + year_str + '.txt'), 'w') as fout:
    #     fout.write('%s: Mean Volume = %0.4f [km3]\n\n' % (which_vol, V/1e9))
    #     for vn in C.keys():
    #         tstr = ' ' + which_vol + ' ' + vn + ' Annual Mean [' + units_dict[vn] + '] '
    #         fout.write(tstr.center(51,'=') + '\n')
    #         for k in Cnorm[vn].keys():
    #             fout.write('%25s %25.3f\n' % (k, Cnorm[vn][k]))
    #
    # pfun.end_plot()

