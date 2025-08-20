"""
Helper script called by make_forcing_main
to generate forcing for point sources 
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import sys
import os
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

#################################################################################
#                   Initialize function and empty dataset                       #
#################################################################################

def make_forcing(N,NT,NRIV,NTRIV,NWWTP_moh,dt_ind, yd_ind,ot_vec,Ldir,enable,trapsP,trapsD,ctag):

    # Start Dataset
    wwtp_ds = xr.Dataset()
    NWWTP = 0
    
    # get year list
    years = [fulldate.year for fulldate in dt_ind]

    # get pre-2020 dates and post-2020 dates:
    dt_ind_thru2020  = dt_ind[dt_ind.year < 2021]
    dt_ind_post2020 = dt_ind[dt_ind.year > 2020]
    yd_ind_thru2020  = [year for date, year in zip(dt_ind, yd_ind) if date.year < 2021]
    yd_ind_post2020 = [year for date, year in zip(dt_ind, yd_ind) if date.year > 2020]

    # adjust date format
    dt_ind = dt_ind.normalize()
    dt_ind_thru2020 = dt_ind_thru2020.normalize()
    dt_ind_post2020 = dt_ind_post2020.normalize()

#################################################################################
#                                Get data                                       #
#################################################################################

    # only get data if WWTPs are enabled
    if enable == True:

        # get raw data
        was24_wwtp_fn = Ldir['data'] / trapsD / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
        was24_wwtp_data_ds = xr.open_dataset(was24_wwtp_fn)

        # first, make sure file exists
        gwi_fn = Ldir['grid'] / 'was24_wwtp_info.csv'
        if not os.path.isfile(gwi_fn):
            print('***Missing was24_wwtp_info.csv file. Please run traps_placement')
            sys.exit()
        # then get the list of WWTPs and indices for this grid
        gwi_df = pd.read_csv(gwi_fn, index_col='rname')
        # if testing, only look at a few sources
        if Ldir['testing']:
            gwi_df = gwi_df.loc[['King County West Point WWTP', 'BELLINGHAM STP'],:]
        # # Test code to simulate WWTP mapped to same grid cell
        # gwi_df = gwi_df.loc[['King County West Point WWTP', 'BELLINGHAM STP'],:]
        # gwi_df.loc['King County West Point WWTP'] = [947.0, 573.0]
        # gwi_df.loc['BELLINGHAM STP'] = [947.0, 573.0]

        # check if there are moh20 and was24 wwtps that are mapped to the same grid cell
        moh20_df = pd.read_csv(Ldir['grid'] / 'moh20_wwtp_info.csv', index_col='rname')
        # # Test code to artifically add overlapping WWTPs
        # gwi_df.loc['FAKE_WAS24_WWTP'] = [947.0, 573.0]
        # moh20_df.loc['FAKE_MOH20_WWTP'] = [947.0, 573.0]
        shared = set(zip(gwi_df.row_py, gwi_df.col_py)) & set(zip(moh20_df.row_py, moh20_df.col_py))
        if shared:
            print("\033[91m\nWARNING: WWTPs in Mohamedali et al. (2020) and Wasielewski et al. (2024) mapped to same grid cell\033[0m")
            # Get matching rows by coordinate
            gwi_matches = gwi_df[gwi_df.set_index(['row_py', 'col_py']).index.isin(shared)]
            moh20_matches = moh20_df[moh20_df.set_index(['row_py', 'col_py']).index.isin(shared)]
            # Print duplicate rnames
            print("Wasielewski et al. (2024) WWTPs:\n", gwi_matches.index.tolist())
            print("Mohamedali et al. (2020) WWTPs:\n", moh20_matches.index.tolist())
            # print next steps
            print('\033[91mPlease write additional code to handle this edge case!!!\033[0m\n')
        else:
            print("WWTPs in Mohamedali et al. (2020) and Wasielewski et al. (2024) mapped to unique grid cells")

        # define directory for point_source climatology
        wwtp_dir = Ldir['LOo'] / 'pre' / trapsP / 'was24_wwtps' /ctag
        traps_type = 'was24wwtp'  

        # get climatological data
        clim_fns = ['Cflow_was24wwtp_fn', 'Ctemp_was24wwtp_fn', 'CDO_was24wwtp_fn',
                    'CNH4_was24wwtp_fn', 'CNO3_was24wwtp_fn', 'CTalk_was24wwtp_fn', 'CTIC_was24wwtp_fn']
        clim_vns = ['flow', 'temp', 'DO', 'NH4', 'NO3', 'Talk', 'TIC']
        for i, clim_fn in enumerate(clim_fns):
            Ldir[clim_fn] = wwtp_dir / 'Data_historical' / ('CLIM_'+clim_vns[i]+'.p')


#################################################################################
#       Combine name of sources that are located at the same grid cell          #
#################################################################################

        # get list of overlapping point sources (i.e. wwtps mapped to same grid cell)
        overlapping_wwtps = gwi_df[gwi_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
        # Remove Lake Stevens 001 and 002 from overlapping list (since they are never active at the same time)
        overlapping_wwtps = [wwtp for wwtp in overlapping_wwtps if wwtp not in ['Lake Stevens 001', 'Lake Stevens 002']]
        # consolidate overlapping point sources
        combined_names = trapsfun.combine_adjacent(overlapping_wwtps)
        gri_df_no_ovrlp = pd.DataFrame(columns=gwi_df.columns) 
        gri_df_no_ovrlp.index.name='rname'
        # loop through original dataframe
        for psname in gwi_df.index: 
            # look for point sources that are in the list of duplicates
            if psname in overlapping_wwtps: 
                # get index in the list of duplicates
                name_index = overlapping_wwtps.index(psname)
                # even index means first occurence of duplicate
                if name_index%2 == 0: 
                    # combine names of duplicates
                    newname = combined_names[int(name_index/2)] 
                    # add combined source to dataframe
                    gri_df_no_ovrlp.loc[newname] = gwi_df.loc[psname] 
                # Note: second duplicate will be dropped
            else:
                # if not a duplicate, then just copy over original info
                gri_df_no_ovrlp.loc[psname] = gwi_df.loc[psname] 

#################################################################################
#                       Prepare dataset for data                                #
#################################################################################

        # get number of wwtps after consolidating overlapping ones
        NWWTP = len(gri_df_no_ovrlp)

        # get the flow, temperature, and nutrient data for these days, is later than 2020
        qtbio_wwtp_df_dict = trapsfun.get_qtbio(gwi_df, dt_ind_post2020, yd_ind_post2020, Ldir, traps_type, trapsD)

        # Add time coordinate
        wwtp_ds['river_time'] = (('river_time',), ot_vec)
        wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        wwtp_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        wwtp_ds['river'] = (('river',), np.arange(NRIV+NTRIV+NWWTP_moh+1,NRIV+NTRIV+NWWTP_moh+NWWTP+1))
        wwtp_ds['river'].attrs['long_name'] = 'marine point source identification number'

        # Add river names
        wwtp_ds['river_name'] = (('river',), list(gri_df_no_ovrlp.index))
        wwtp_ds['river_name'].attrs['long_name'] = 'point source name'

#################################################################################
#  Add vertical distribution of sources. All WWTPs discharge from bottom layer  #
#################################################################################

        # Add Vshape
        vn = 'river_Vshape'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = ('s_rho', 'river')
        # All discharge coming from the bottom layer
        Vshape = np.zeros((N, NWWTP))
        Vshape[0,:] = 1
        wwtp_ds[vn] = (dims, Vshape)
        wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

#################################################################################
#             Add indices of sources. WWTPs located on the rho-grid             #
#################################################################################

        # Add position and direction
        for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            # set point source diretion to enter vertically (Dsrc = 2)
            if vn == 'river_direction':
                wwtp_direction = 2 * np.ones(NWWTP) 
                wwtp_ds[vn] = (('river',), wwtp_direction)
            # Add X-position (column index)
            elif vn == 'river_Xposition':
                X_vec = np.nan * np.ones(NWWTP)
                for ii,wn in enumerate(gri_df_no_ovrlp.index):
                    X_vec[ii] = gri_df_no_ovrlp.loc[wn, 'col_py']
                wwtp_ds[vn] = (('river',), X_vec)
            # Add E-position (row index)
            elif vn == 'river_Eposition':
                E_vec = np.nan * np.ones(NWWTP)
                for ii,wn in enumerate(gri_df_no_ovrlp.index):
                    E_vec[ii] = gri_df_no_ovrlp.loc[wn, 'row_py']
                wwtp_ds[vn] = (('river',), E_vec)
            # add metadata
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#                               Add source flowrate                             #
#################################################################################

        # Add transport
        vn = 'river_transport'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('river',)
        Q_mat = np.zeros((NT, NWWTP))
        for rr,rn in enumerate(gri_df_no_ovrlp.index):
            # consolidate sources located at same grid cell (sum flowrates)
            if '+' in rn:
                # split into individual point sources
                [wwtp1,wwtp2] = rn.split('+')
                # get individual point source flowrates

                # get flow through 2020 from actual data
                if len(dt_ind_thru2020) == 0:
                    flow1_thru2020 = np.array([])
                else:
                    flow1_thru2020 = was24_wwtp_data_ds.flow.sel(source=
                        was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == wwtp1].item(),
                        date=dt_ind_thru2020).values     
                # get flow post 2020 from climatologies         
                flow1_post2020 = qtbio_wwtp_df_dict[wwtp1]['flow'].values
                # concatenate the flows
                flow1 = np.concatenate([flow1_thru2020, flow1_post2020])

                # get flow through 2020 from actual data
                if len(dt_ind_thru2020) == 0:
                    flow2_thru2020 = np.array([])
                else:
                    flow2_thru2020 = was24_wwtp_data_ds.flow.sel(source=
                        was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == wwtp2].item(),
                        date=dt_ind_thru2020).values
                # get flow post 2020 from climatologies         
                flow2_post2020 = qtbio_wwtp_df_dict[wwtp2]['flow'].values
                # concatenate the flows
                flow2 = np.concatenate([flow2_thru2020, flow2_post2020])
                
                # combine point source flow
                flow = flow1 + flow2
            else:
                # get flow through 2020 from actual data
                if len(dt_ind_thru2020) == 0:
                    flow_thru2020 = np.array([])
                else:
                    flow_thru2020 = was24_wwtp_data_ds.flow.sel(source=
                        was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                        date=dt_ind_thru2020).values
                # get flow post 2020 from climatologies
                flow_post2020 = qtbio_wwtp_df_dict[rn]['flow'].values
                # concatenate the flows
                flow = np.concatenate([flow_thru2020, flow_post2020])
                
            # update flowrate
            Q_mat[:,rr] = flow
        # add metadata
        wwtp_ds[vn] = (dims, Q_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                         Add source salinity and temp                          #
#################################################################################

        # Add salinity and temperature
        for vn in ['river_salt', 'river_temp']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            # salinity is always zero
            if vn == 'river_salt':
                TS_mat = np.zeros((NT, N, NWWTP))
            # get temperature from climatology
            elif vn == 'river_temp':
                TS_mat = np.nan * np.zeros((NT, N, NWWTP))
                for rr,rn in enumerate(gri_df_no_ovrlp.index):
                    # consolidate sources located at same grid cell
                    if '+' in rn:
                        # split into individual point sources
                        [wwtp1,wwtp2] = rn.split('+')
                        # calculate weighted average (based on flowrate)
                        # get temp through 2020 from actual data
                        if len(dt_ind_thru2020) == 0:
                            temps_thru2020 = np.array([])
                        else:
                            temps_thru2020 = trapsfun.weighted_average_ds('temp', dt_ind_thru2020, was24_wwtp_data_ds, wwtp1, wwtp2)
                        # get temp post 2020 from climatologies
                        temps_post2020 = trapsfun.weighted_average('temp',qtbio_wwtp_df_dict[wwtp1], qtbio_wwtp_df_dict[wwtp2])
                        # concatenate the temps
                        temps = np.concatenate([temps_thru2020, temps_post2020])

                    else:
                        # get temp through 2020 from actual data
                        if len(dt_ind_thru2020) == 0:
                            temps_thru2020 = np.array([])
                        else:
                            temps_thru2020 = was24_wwtp_data_ds.temp.sel(source=
                                was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                                date=dt_ind_thru2020).values
                        # get temp post 2020 from climatologies
                        temps_post2020 = qtbio_wwtp_df_dict[rn]['temp'].values
                        # concatenate the temps
                        temps = np.concatenate([temps_thru2020, temps_post2020])
                        
                    print('-- {}: filled from raw Wasielewski et al. (2024) dataset'.format(rn))
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in was24 wwtp river_temp!')
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, TS_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                            Add source biology                                 #
#################################################################################

        # Add biologeochemistry parameters
        for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
        # for var in ['NO3', 'NH4', 'TIC', 'Talk', 'DO']:
            var_post = var
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            for rr,rn in enumerate(gri_df_no_ovrlp.index):
                # adjust names to get data from dataset
                if var == 'TAlk':
                    var = 'Talk'
                if var == 'Oxyg':
                    var = 'DO'
                # consolidate sources located at same grid cell
                if '+' in rn:
                    # split into individual point sources
                    [wwtp1,wwtp2] = rn.split('+')
                    # calculate weighted average (based on flowrate)
                    # get values through 2020 from actual data
                    if len(dt_ind_thru2020) == 0:
                        bvals_thru2020 = np.array([])
                    else:
                        bvals_thru2020 = trapsfun.weighted_average_ds(var, dt_ind_thru2020, was24_wwtp_data_ds, wwtp1, wwtp2)
                    # get values post 2020 from climatologies
                    bvals_post2020 = trapsfun.weighted_average(var_post,qtbio_wwtp_df_dict[wwtp1], qtbio_wwtp_df_dict[wwtp2])
                    # concatenate the values
                    bvals = np.concatenate([bvals_thru2020, bvals_post2020])
                else:
                    # get values through 2020 from actual data
                    if len(dt_ind_thru2020) == 0:
                        bvals_thru2020 = np.array([])
                    else:
                        bvals_thru2020 = was24_wwtp_data_ds[var].sel(
                                source=was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                                date=dt_ind_thru2020)
                    # get values post 2020 from climatologies
                    bvals_post2020 = qtbio_wwtp_df_dict[rn][var_post].values
                    # concatenate the values
                    bvals = np.concatenate([bvals_thru2020, bvals_post2020])
                    
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
            # check for nans
            if np.isnan(B_mat).any():
                print('Error from traps: nans in was24 wwtp bio!')
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                  All other biology variables are zero                         #
#################################################################################

        # Add remaining biology (see the lineup near the end of fennel_var.h)
        # Right now, this is simply filling everything with zeros
        bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
        for bvn in bvn_list:
            vn = 'river_' + bvn
            vinfo = zrfun.get_varinfo(vn)
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            # loop through all sources and fill with zeros
            for rr,rn in enumerate(gri_df_no_ovrlp.index):
                for nn in range(N):
                    B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
            # check for nans
            if np.isnan(B_mat).any():
                print('Error from traps: nans in B_mat for was24 wwtp ' + vn)
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#          Return WWTP forcing dataset in the form that ROMS expects            #
#################################################################################

    # print(wwtp_ds.river_transport)
    # print(wwtp_ds.river_NH4)

    return wwtp_ds, NWWTP
