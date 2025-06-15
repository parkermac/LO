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

def make_forcing(N,NT,NRIV,NTRIV,dt_ind, yd_ind,ot_vec,Ldir,enable,trapsP,trapsD,ctag):

    # Start Dataset
    wwtp_ds = xr.Dataset()
    NWWTP = 0

    # get year list
    years = [fulldate.year for fulldate in dt_ind]

#################################################################################
#                                Get data                                       #
#################################################################################

    # only get data if WWTPs are enabled
    if enable == True:

        # define directory for point_source climatology
        wwtp_dir = Ldir['LOo'] / 'pre' / trapsP / 'moh20_wwtps' /ctag
        traps_type = 'wwtp'  

        # get climatological data
        clim_fns = ['Cflow_wwtp_fn', 'Ctemp_wwtp_fn', 'CDO_wwtp_fn',
                    'CNH4_wwtp_fn', 'CNO3_wwtp_fn', 'CTalk_wwtp_fn', 'CTIC_wwtp_fn']
        clim_vns = ['flow', 'temp', 'DO', 'NH4', 'NO3', 'Talk', 'TIC']
        for i, clim_fn in enumerate(clim_fns):
            Ldir[clim_fn] = wwtp_dir / 'Data_historical' / ('CLIM_'+clim_vns[i]+'.p')

        # first, make sure file exists
        gwi_fn = Ldir['grid'] / 'moh20_wwtp_info.csv'
        if not os.path.isfile(gwi_fn):
            print('***Missing moh20_wwtp_info.csv file. Please run traps_placement')
            sys.exit()
        # then get the list of point sources and indices for this grid
        gwi_df = pd.read_csv(gwi_fn, index_col='rname')
        # if testing, only look at a few sources
        if Ldir['testing']:
            gwi_df = gwi_df.loc[['LOTT', 'Iona'],:]
        
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

        # get the flow, temperature, and nutrient data for these days
        qtbio_wwtp_df_dict = trapsfun.get_qtbio(gwi_df, dt_ind, yd_ind, Ldir, traps_type, trapsD)

        # Add time coordinate
        wwtp_ds['river_time'] = (('river_time',), ot_vec)
        wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        wwtp_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        wwtp_ds['river'] = (('river',), np.arange(NRIV+NTRIV+1,NRIV+NTRIV+NWWTP+1))
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
                qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                flow1 = qtbio_wwtp_df_1['flow'].values
                flow2 = qtbio_wwtp_df_2['flow'].values
                # combine point source flow
                flow = flow1 + flow2
            else:
                qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                flow = qtbio_wwtp_df['flow'].values
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
                        # get individual point source dataframe
                        qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                        qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                        # calculate weighted average (based on flowrate)
                        temps = trapsfun.weighted_average('temp',qtbio_wwtp_df_1, qtbio_wwtp_df_2)
                    else:
                        qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                        temps = qtbio_wwtp_df['temp'].values
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
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            for rr,rn in enumerate(gri_df_no_ovrlp.index):
                # consolidate sources located at same grid cell
                if '+' in rn:
                    # split into individual point sources
                    [wwtp1,wwtp2] = rn.split('+')
                    # get individual point source dataframe
                    qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                    qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                    # calculate weighted average (based on flowrate)
                    bvals = trapsfun.weighted_average(var,qtbio_wwtp_df_1, qtbio_wwtp_df_2)
                else:
                    qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                    bvals = qtbio_wwtp_df[var].values
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

    return wwtp_ds, NWWTP
