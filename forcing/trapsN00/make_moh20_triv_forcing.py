"""
Helper script called by make_forcing_main
to generate forcing for tiny rivers 
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

def make_forcing(N,NT,NRIV,dt_ind, yd_ind,ot_vec,Ldir,enable,trapsP,trapsD,ctag):
    # Start Dataset
    triv_ds = xr.Dataset()
    NTRIV = 0

#################################################################################
#                                Get data                                       #
#################################################################################

    # only get data if tiny rivers are enabled
    if enable == True:

        # define directory for tiny river climatology
        triv_dir = Ldir['LOo'] / 'pre' / trapsP / 'moh20_tinyrivers' / ctag
        traps_type = 'triv'

        # get climatological data
        clim_fns = ['Cflow_triv_fn', 'Ctemp_triv_fn', 'CDO_triv_fn',
                    'CNH4_triv_fn', 'CNO3_triv_fn', 'CTalk_triv_fn', 'CTIC_triv_fn']
        clim_vns = ['flow', 'temp', 'DO', 'NH4', 'NO3', 'Talk', 'TIC']
        for i, clim_fn in enumerate(clim_fns):
            Ldir[clim_fn] = triv_dir / 'Data_historical' / ('CLIM_'+clim_vns[i]+'.p')

        # first, make sure file exists
        gtri_fn = Ldir['grid'] / 'moh20_triv_info.csv'
        if not os.path.isfile(gtri_fn):
            print('***Missing moh20_triv_info.csv file. Please run traps_placement')
            sys.exit()
        # then get the list of tiny rivers and indices for this grid
        gtri_df = pd.read_csv(gtri_fn, index_col='rname')
        # if testing, only look at a few sources
        if Ldir['testing']:
            gtri_df = gtri_df.loc[['Birch Bay', 'Purdy Cr', 'Burley Cr', 'Perry Cr','McLane Cr'],:]
        
#################################################################################
#       Combine name of sources that are located at the same grid cell          #
#################################################################################

        # get list of overlapping rivers (i.e. rivers mapped to same grid cell)
        overlapping_trivs = gtri_df[gtri_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
        # consolidate overlapping rivers
        combined_names = trapsfun.combine_adjacent(overlapping_trivs)
        gri_df_no_ovrlp = pd.DataFrame(columns=gtri_df.columns) 
        gri_df_no_ovrlp.index.name='rname'
        # loop through original dataframe
        for trname in gtri_df.index:
            # look for rivers that are in the list of duplicates
            if trname in overlapping_trivs: 
                # get index in the list of duplicates
                name_index = np.where(overlapping_trivs == trname)[0][0]
                # even index means first occurence of duplicate
                if name_index%2 == 0: 
                    # combine names of duplicates
                    newname = combined_names[int(name_index/2)] 
                    # add combined source to dataframe
                    gri_df_no_ovrlp.loc[newname] = gtri_df.loc[trname]
                # Note: second duplicate will be dropped
                # (so idir, isign, and uv will come from the first duplicate)
            else:
                # if not a duplicate, then just copy over original info
                gri_df_no_ovrlp.loc[trname] = gtri_df.loc[trname]

#################################################################################
#                       Prepare dataset for data                                #
#################################################################################

        # get number of tiny rivers after consolidating overlapping ones
        NTRIV = len(gri_df_no_ovrlp)

        # get the flow, temperature, and nutrient data for these days
        qtbio_triv_df_dict = trapsfun.get_qtbio(gtri_df, dt_ind, yd_ind, Ldir, traps_type, trapsD)

        # Add time coordinate
        triv_ds['river_time'] = (('river_time',), ot_vec)
        triv_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        triv_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        triv_ds['river'] = (('river',), np.arange(NRIV+1,NRIV+NTRIV+1))
        triv_ds['river'].attrs['long_name'] = 'tiny river runoff identification number'

#################################################################################
#  Add vertical distribution of sources. All rivers discharge uniformly in z    #
#################################################################################

        # Add Vshape
        vn = 'river_Vshape'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = ('s_rho', 'river')
        # For Vtransform = 2, even spacing is a good approximation, and
        # we implement this by using 1/N as the fraction in each vertical cell.
        Vshape = (1/N) * np.ones((N, NTRIV))
        triv_ds[vn] = (dims, Vshape)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#            Add indices of sources. Rivers located on u- or v-grid             #
#################################################################################

        # Add position and direction
        for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            # get river direction (idir)
            if vn == 'river_direction':
                triv_ds[vn] = (('river',), gri_df_no_ovrlp.idir.to_numpy().astype(int))
            # Add X-position (column index on u-grid)
            elif vn == 'river_Xposition':
                X_vec = np.nan * np.ones(NTRIV)
                for ii,rn in enumerate(gri_df_no_ovrlp.index):
                    if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                        X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py'] + 1
                    elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                        X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py']
                triv_ds[vn] = (('river',), X_vec)
            # Add E-position (row index on v-grid)
            elif vn == 'river_Eposition':
                E_vec = np.nan * np.ones(NTRIV)
                for ii,rn in enumerate(gri_df_no_ovrlp.index):
                    if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                        E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py']
                    elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                        E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py'] + 1
                triv_ds[vn] = (('river',), E_vec)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#                               Add source flowrate                             #
#################################################################################

        # Add transport
        vn = 'river_transport'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('river',)
        Q_mat = np.zeros((NT, NTRIV))
        for rr,rn in enumerate(gri_df_no_ovrlp.index):
            # consolidate sources located at same grid cell (sum flowrates)
            if '+' in rn:
                # split into individual rivers
                [riv1,riv2] = rn.split('+')
                # get individual river flowrates
                qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                flow1 = qtbio_triv_df_1['flow'].values
                flow2 = qtbio_triv_df_2['flow'].values
                # combine river flow
                flow = flow1 + flow2
            else:
                qtbio_triv_df = qtbio_triv_df_dict[rn]
                flow = qtbio_triv_df['flow'].values
            Q_mat[:,rr] = flow * gri_df_no_ovrlp.loc[rn, 'isign']
        # add metadata
        triv_ds[vn] = (dims, Q_mat)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']
        triv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                         Add source salinity and temp                          #
#################################################################################

        # Add salinity and temperature
        for vn in ['river_salt', 'river_temp']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            # salinity is always zero
            if vn == 'river_salt':
                TS_mat = np.zeros((NT, N, NTRIV))
            # get temperature from climatology
            elif vn == 'river_temp':
                TS_mat = np.nan * np.zeros((NT, N, NTRIV))
                for rr,rn in enumerate(gri_df_no_ovrlp.index):
                    # consolidate sources located at same grid cell
                    if '+' in rn:
                        # split into individual rivers
                        [riv1,riv2] = rn.split('+')
                        # get individual river dataframe
                        qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                        qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                        # calculate weighted average (based on flowrate)
                        temps = trapsfun.weighted_average('temp',qtbio_triv_df_1, qtbio_triv_df_2)
                    else:
                        qtbio_triv_df = qtbio_triv_df_dict[rn]
                        temps = qtbio_triv_df['temp'].values
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river river_temp!')
                sys.exit()
            # add metadata
            triv_ds[vn] = (dims, TS_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                            Add source biology                                 #
#################################################################################

        # Add biologeochemistry parameters
        for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NTRIV))
            for rr,rn in enumerate(gri_df_no_ovrlp.index):
                # consolidate sources located at same grid cell
                if '+' in rn:
                    # split into individual rivers
                    [riv1,riv2] = rn.split('+')
                    # get individual river dataframe
                    qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                    qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                    # calculate weighted average (based on flowrate)
                    bvals = trapsfun.weighted_average(var,qtbio_triv_df_1, qtbio_triv_df_2)
                else:
                    qtbio_triv_df = qtbio_triv_df_dict[rn]
                    bvals = qtbio_triv_df[var].values
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
            # check for nans
            if np.isnan(B_mat).any():
                print('Error from traps: nans in tiny river bio!')
                sys.exit()
            # add metadata
            triv_ds[vn] = (dims, B_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

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
            B_mat = np.nan * np.zeros((NT, N, NTRIV))
            # loop through all sources and fill with zeros
            for rr,rn in enumerate(gri_df_no_ovrlp.index):
                for nn in range(N):
                    B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
            # check for nans
            if np.isnan(B_mat).any():
                print('Error from traps: nans in B_mat for tiny river ' + vn)
                sys.exit()
            # add metadata
            triv_ds[vn] = (dims, B_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                                     Add river names                           #
#################################################################################

        # Rename rivers that share name with WWTP. This code appends ' R' at the end of the river name
        duplicates = ['Port Angeles', 'Port Townsend', 'Birch Bay', 'Port Gamble', 'Gig Harbor']
        gri_df_no_ovrlp.index = np.where(gri_df_no_ovrlp.index.isin(duplicates),
                                         gri_df_no_ovrlp.index + ' R', gri_df_no_ovrlp.index)

        # Add river names
        triv_ds['river_name'] = (('river',), list(gri_df_no_ovrlp.index))
        triv_ds['river_name'].attrs['long_name'] = 'tiny river name'

#################################################################################
#          Return triv forcing dataset in the form that ROMS expects            #
#################################################################################

    return triv_ds, NTRIV