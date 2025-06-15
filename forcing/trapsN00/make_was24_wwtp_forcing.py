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

def make_forcing(N,NT,NRIV,NTRIV,NWWTP_moh,dt_ind, yd_ind,ot_vec,Ldir,enable,trapsD):

    # Start Dataset
    wwtp_ds = xr.Dataset()
    NWWTP = 0

    # get year list
    years = [fulldate.year for fulldate in dt_ind]

    # adjust date format
    dt_ind = dt_ind.normalize()

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
                flow1 = was24_wwtp_data_ds.flow.sel(source=
                    was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == wwtp1].item(),
                    date=dt_ind).values
                flow2 = was24_wwtp_data_ds.flow.sel(source=
                    was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == wwtp2].item(),
                    date=dt_ind).values
                # combine point source flow
                flow = flow1 + flow2
            else:
                flow = was24_wwtp_data_ds.flow.sel(source=
                    was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                    date=dt_ind).values
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
                        temps = trapsfun.weighted_average_ds('temp', dt_ind, was24_wwtp_data_ds, wwtp1, wwtp2)
                    else:
                        temps = was24_wwtp_data_ds.temp.sel(source=
                            was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                            date=dt_ind).values
                    print('-- {}: filled from raw Wasielewski et al. (2024) dataset'.format(rn))
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in point source river_temp!')
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
                    bvals = trapsfun.weighted_average_ds(var, dt_ind, was24_wwtp_data_ds, wwtp1, wwtp2)
                else:
                    bvals = was24_wwtp_data_ds[var].sel(
                            source=was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                            date=dt_ind)
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river bio!')
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
                print('Error from traps: nans in B_mat for tiny river ' + vn)
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
