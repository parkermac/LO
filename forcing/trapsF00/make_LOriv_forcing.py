"""
Helper script called by make_forcing_main
to generate forcing for pre-existing LO rivers 
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

def make_forcing(N,NT,dt_ind,yd_ind,ot_vec,dt1,days,Ldir,trapsP,trapsD):
    # Start Dataset
    LOriv_ds = xr.Dataset()

#################################################################################
#                                Get data                                       #
#################################################################################

    # Load a dataframe with info for rivers to get
    if Ldir['gridname'] == 'cas7':
        ctag = 'lo_base'
    else:
        print('You need to specify a gridname for this ctag.')
        sys.exit()

    # define directory for pre-existing LO river climatology
    ri_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag
    ri_df_fn = ri_dir / 'river_info.p'
    ri_df = pd.read_pickle(ri_df_fn)

    # get historical and climatological data files
    Ldir['Hflow_fn'] = ri_dir / 'Data_historical' / ('ALL_flow.p')
    Ldir['Cflow_fn'] = ri_dir / 'Data_historical' / ('CLIM_flow.p')
    Ldir['Ctemp_fn'] = ri_dir / 'Data_historical' / ('CLIM_temp.p')

    # define directory for pre-existing LO river bio climatology
    LObio_dir = Ldir['LOo'] / 'pre' / trapsP / 'LO_rivbio' / ctag
    traps_type = 'LOriv'

    # get climatological data
    clim_fns = ['CDO_LOriv_fn', 'CNH4_LOriv_fn', 'CNO3_LOriv_fn',
                'CTalk_LOriv_fn', 'CTIC_LOriv_fn']
    clim_vns = ['DO', 'NH4', 'NO3', 'Talk', 'TIC']
    for i, clim_fn in enumerate(clim_fns):
        Ldir[clim_fn] = LObio_dir / 'Data_historical' / ('CLIM_'+clim_vns[i]+'.p')

    # get names of rivers for which Ecology has biogeochem data
    # these are the names the LiveOcean calls them.
    # Later, they will be converted to the name Ecology/SSM uses
    repeatrivs_fn = Ldir['data'] / trapsD / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    LObio_names_all = list(repeatrivs_df.loc[repeatrivs_df['in_both'] == 1, 'LO_rname'])
    # remove the weird rivers
    weird_rivers = ['Alberni Inlet', 'Chehalis R', 'Gold River',
                    'Willapa R', 'Columbia R', 'Comox']
    # These are the names that LO uses
    LObio_names = [rname for rname in LObio_names_all if trapsfun.LO2SSM_name(rname,trapsD) not in weird_rivers]

    # get the list of rivers and indices for this grid
    gri_fn = Ldir['grid'] / 'river_info.csv'
    gri_df = pd.read_csv(gri_fn, index_col='rname')
    if Ldir['testing']:
        gri_df = gri_df.loc[['columbia', 'skagit'],:]

#################################################################################
#                       Prepare dataset for data                                #
#################################################################################

    # get number of pre-existing LO rivers
    NRIV = len(gri_df)

    # associate rivers with ones that have temperature climatology data
    ri_df = rivfun.get_tc_rn(ri_df)

    # get the flow and temperature data for these days
    qt_df_dict = rivfun.get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days)
    # get the biology for LO pre-existing rivers for which Ecology has data
    LObio_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type, trapsD)

    # Add time coordinate
    LOriv_ds['river_time'] = (('river_time',), ot_vec)
    LOriv_ds['river_time'].attrs['units'] = Lfun.roms_time_units
    LOriv_ds['river_time'].attrs['long_name'] = 'river time'

    # Add river coordinate
    LOriv_ds['river'] = (('river',), np.arange(1,NRIV+1))
    LOriv_ds['river'].attrs['long_name'] = 'river runoff identification number'

    # Add river names
    LOriv_ds['river_name'] = (('river',), list(gri_df.index))
    LOriv_ds['river_name'].attrs['long_name'] = 'river name'

#################################################################################
#  Add vertical distribution of sources. All rivers discharge uniformly in z    #
#################################################################################

    # Add Vshape
    vn = 'river_Vshape'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = ('s_rho', 'river')
    # For Vtransform = 2, even spacing is a good approximation, and
    # we implement this by using 1/N as the fraction in each vertical cell.
    Vshape = (1/N) * np.ones((N, NRIV))
    LOriv_ds[vn] = (dims, Vshape)
    LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#            Add indices of sources. Rivers located on u- or v-grid             #
#################################################################################

    # Add position and direction
    for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        # get river direction (idir)
        if vn == 'river_direction':
            LOriv_ds[vn] = (('river',), gri_df.idir.to_numpy())
        # Add X-position (column index on v-grid)
        elif vn == 'river_Xposition':
            X_vec = np.nan * np.ones(NRIV)
            for ii,rn in enumerate(gri_df.index):
                if gri_df.loc[rn, 'idir'] == 0:
                    X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
                elif gri_df.loc[rn, 'idir'] == 1:
                    X_vec[ii] = gri_df.loc[rn, 'col_py']
            LOriv_ds[vn] = (('river',), X_vec)
        # Add E-position (row index on u-grid)
        elif vn == 'river_Eposition':
            E_vec = np.nan * np.ones(NRIV)
            for ii,rn in enumerate(gri_df.index):
                if gri_df.loc[rn, 'idir'] == 0:
                    E_vec[ii] = gri_df.loc[rn, 'row_py']
                elif gri_df.loc[rn, 'idir'] == 1:
                    E_vec[ii] = gri_df.loc[rn, 'row_py'] + 1
            LOriv_ds[vn] = (('river',), E_vec)
        LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#                               Add source flowrate                             #
#################################################################################

    # Add transport
    vn = 'river_transport'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('river',)
    Q_mat = np.zeros((NT, NRIV))
    for rr,rn in enumerate(gri_df.index):
        qt_df = qt_df_dict[rn]
        flow = qt_df['final'].values
        Q_mat[:,rr] = flow * gri_df.loc[rn, 'isign']
    # add metadata
    LOriv_ds[vn] = (dims, Q_mat)
    LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
    LOriv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                         Add source salinity and temp                          #
#################################################################################

    # Add salinity and temperature
    for vn in ['river_salt', 'river_temp']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('s_rho', 'river')
        # salinity is always zero
        if vn == 'river_salt':
            TS_mat = np.zeros((NT, N, NRIV))
        # get temperature
        elif vn == 'river_temp':
            TS_mat = np.nan * np.zeros((NT, N, NRIV))
            for rr,rn in enumerate(gri_df.index):
                qt_df = qt_df_dict[rn]
                for nn in range(N):
                    TS_mat[:, nn, rr] = qt_df['temperature'].values
        # check for nans
        if np.isnan(TS_mat).any():
            print('Error from pre-existing LO rivers: nans in river_temp!')
            sys.exit()
        # add metadata
        LOriv_ds[vn] = (dims, TS_mat)
        LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
        LOriv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                            Add source biology                                 #
#################################################################################
   
    # Add biology (see the lineup near the end of fennel_var.h)
    bvn_list = ['NO3', 'NH4', 'Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo',
            'TIC', 'TAlk', 'LDeC', 'SDeC', 'Oxyg']
    for bvn in bvn_list:
        vn = 'river_' + bvn
        vinfo = zrfun.get_varinfo(vn)
        dims = (vinfo['time'],) + ('s_rho', 'river')
        B_mat = np.nan * np.zeros((NT, N, NRIV))
        for rr,rn in enumerate(gri_df.index):
            # Add biogeochem climatology for rivers for which Ecology have data
            if rn in LObio_names and bvn in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
                # get corresponding Ecology/SSM river name
                rn_SSM = trapsfun.LO2SSM_name(rn, trapsD)
                # get the biogeochem values from climatology
                bio_LOriv_df = LObio_df_dict[rn_SSM]
                bvals = bio_LOriv_df[bvn].values
                # use ammonium values from Susan Allen for Fraser
                # https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?
                # doi=10.1029%2F2019JC015766&file=jgrc24099-sup-0001-Text_SI-S01.pdf
                if rn == 'fraser' and bvn == 'NH4':
                    bvals = 4.43 * np.ones(NT) # uM = mmol/m3
            # If Ecology doesn't have data, use default LO bio
            else:
                bvals = rivfun.get_bio_vec(bvn, rn, yd_ind)
            for nn in range(N):
                B_mat[:, nn, rr] = bvals
        # check for nans
        if np.isnan(B_mat).any():
            print('Error from pre-existing LO rivers: nans in B_mat for ' + vn)
            sys.exit()
        # add metadata
        LOriv_ds[vn] = (dims, B_mat)
        LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
        LOriv_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#          Return LOriv forcing dataset in the form that ROMS expects           #
#################################################################################

    return LOriv_ds, NRIV