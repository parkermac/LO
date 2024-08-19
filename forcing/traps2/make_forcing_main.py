"""
This is the main program for making the RIVER and TRAPS forcing file, for the
updated ROMS

Test on pc in ipython:
run make_forcing_main.py -g cas6 -r backfill -d 2021.01.01 -f traps2 -test True

2023.04.11 Updated to work with pre/river1 data format.

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ENABLE OR DISABLE TINY RIVERS AND/OR POINT SOURCES
enable_tinyrivers = True 
enable_pointsources = True

# ****************** CASE-SPECIFIC CODE *****************

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + date_string) / Ldir['frc']

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)
    reload(rivfun)

out_fn = out_dir / 'rivers.nc'
out_fn.unlink(missing_ok=True)

# set up the time index for the record
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
yd_ind = pd.Index(dt_ind.dayofyear)

ot_vec = np.array([Lfun.datetime_to_modtime(item) for item in dt_ind])
NT = len(ot_vec)

S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
N = S['N']

grid_fn = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(grid_fn, only_G=True)

###############################################################################################
# LIVEOCEAN PRE-EXISTING RIVERS

# Load a dataframe with info for rivers to get
if Ldir['gridname'] == 'cas6':
    ctag = 'lo_base'
else:
    print('You need to specify a gridname for this ctag.')
    sys.exit()

ri_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag
ri_df_fn = ri_dir / 'river_info.p'
ri_df = pd.read_pickle(ri_df_fn)

# get historical and climatological data files
Ldir['Hflow_fn'] = ri_dir / 'Data_historical' / ('ALL_flow.p')
Ldir['Cflow_fn'] = ri_dir / 'Data_historical' / ('CLIM_flow.p')
Ldir['Ctemp_fn'] = ri_dir / 'Data_historical' / ('CLIM_temp.p')

# get biologeochem data for rivers for which Ecology has data
LObio_dir = Ldir['LOo'] / 'pre' / 'traps' / 'LO_rivbio'
traps_type = 'LOriv'
# climatological data files
year0 = 1999
year1 = 2017
# climatological data
Ldir['CDO_LOriv_fn']   = LObio_dir / 'Data_historical' / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p')
Ldir['CNH4_LOriv_fn']  = LObio_dir / 'Data_historical' / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p')
Ldir['CNO3_LOriv_fn']  = LObio_dir / 'Data_historical' / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p')
Ldir['CTalk_LOriv_fn'] = LObio_dir / 'Data_historical' / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p')
Ldir['CTIC_LOriv_fn']  = LObio_dir / 'Data_historical' / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p')

# get names of rivers for which Ecology has biogeochem data
# these are the names the LiveOcean calls them. Later, they will be converted to the name Ecology/SSM uses
repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
repeatrivs_df = pd.read_excel(repeatrivs_fn)
LObio_names_all = list(repeatrivs_df.loc[repeatrivs_df['in_both'] == 1, 'LO_rname'])
# remove the weird rivers
weird_rivers = ['Alberni Inlet', 'Chehalis R', 'Gold River', 'Willapa R', 'Columbia R', 'Comox']
# These are the names that LO uses
LObio_names = [rname for rname in LObio_names_all if trapsfun.LO2SSM_name(rname) not in weird_rivers]

# get the list of rivers and indices for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')
if Ldir['testing']:
    gri_df = gri_df.loc[['columbia', 'skagit'],:]
NRIV = len(gri_df)

# associate rivers with ones that have temperature climatology data
ri_df = rivfun.get_tc_rn(ri_df)

# get the flow and temperature data for these days
qt_df_dict = rivfun.get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days)
# get the biology for LO pre-existing rivers for which Ecology has data
LObio_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type)

# Start Dataset
LOriv_ds = xr.Dataset()

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

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NRIV))
LOriv_ds[vn] = (dims, Vshape)
LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        LOriv_ds[vn] = (('river',), gri_df.idir.to_numpy().astype(int))
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NRIV)
        ii = 0
        for rn in gri_df.index:
            if gri_df.loc[rn, 'idir'] == 0:
                X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
            elif gri_df.loc[rn, 'idir'] == 1:
                X_vec[ii] = gri_df.loc[rn, 'col_py']
            ii += 1
        LOriv_ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NRIV)
        ii = 0
        for rn in gri_df.index:
            if gri_df.loc[rn, 'idir'] == 0:
                E_vec[ii] = gri_df.loc[rn, 'row_py']
            elif gri_df.loc[rn, 'idir'] == 1:
                E_vec[ii] = gri_df.loc[rn, 'row_py'] + 1
            ii += 1
        LOriv_ds[vn] = (('river',), E_vec)
    LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
        

# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NRIV))
rr = 0
for rn in gri_df.index:
    qt_df = qt_df_dict[rn]
    flow = qt_df['final'].values
    Q_mat[:,rr] = flow * gri_df.loc[rn, 'isign']
    rr += 1
LOriv_ds[vn] = (dims, Q_mat)
LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
LOriv_ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature
for vn in ['river_salt', 'river_temp']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')
    if vn == 'river_salt':
        TS_mat = np.zeros((NT, N, NRIV))
    elif vn == 'river_temp':
        TS_mat = np.nan * np.zeros((NT, N, NRIV))
        rr = 0
        for rn in gri_df.index:
            qt_df = qt_df_dict[rn]
            for nn in range(N):
                TS_mat[:, nn, rr] = qt_df['temperature'].values
            rr += 1
    if np.isnan(TS_mat).any():
        print('Error from riv00: nans in river_temp!')
        sys.exit()
    LOriv_ds[vn] = (dims, TS_mat)
    LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
    LOriv_ds[vn].attrs['units'] = vinfo['units']
    
# Add biology (see the lineup near the end of fennel_var.h)
bvn_list = ['NO3', 'NH4', 'Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo',
        'TIC', 'TAlk', 'LDeC', 'SDeC', 'Oxyg']
for bvn in bvn_list:
    vn = 'river_' + bvn
    vinfo = zrfun.get_varinfo(vn)
    dims = (vinfo['time'],) + ('s_rho', 'river')
    B_mat = np.nan * np.zeros((NT, N, NRIV))
    rr = 0
    for rn in gri_df.index:
        # Add biogeochem climatology for rivers for which Ecology have data
        if rn in LObio_names and bvn in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
            # get corresponding Ecology/SSM river name
            rn_SSM = trapsfun.LO2SSM_name(rn)
            # get the biogeochem values from climatology
            bio_LOriv_df = LObio_df_dict[rn_SSM]
            bvals = bio_LOriv_df[bvn].values
        # If Ecology doesn't have data, use default LO bio
        else:
            bvals = rivfun.get_bio_vec(bvn, rn, yd_ind)
        for nn in range(N):
            B_mat[:, nn, rr] = bvals
        rr += 1
    if np.isnan(B_mat).any():
        print('Error from riv00: nans in B_mat for ' + vn)
        sys.exit()
    LOriv_ds[vn] = (dims, B_mat)
    LOriv_ds[vn].attrs['long_name'] = vinfo['long_name']
    LOriv_ds[vn].attrs['units'] = vinfo['units']
###########################################################################################
# FORCING FOR TINY RIVERS

# Start Dataset
triv_ds = xr.Dataset()

NTRIV = 0

if enable_tinyrivers == True:

    # Run placement algorithm to put tiny rivers on LiveOcean grid
    trapsfun.traps_placement('riv')

    # define directory for tiny river climatology
    tri_dir = Ldir['LOo'] / 'pre' / 'traps' / 'tiny_rivers'
    traps_type = 'triv'

    # climatological data files
    year0 = 1999
    year1 = 2017
    # climatological data
    Ldir['Cflow_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['Ctemp_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CDO_triv_fn']   = tri_dir / 'Data_historical' / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CNH4_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CNO3_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CTalk_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CTIC_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p')

    # get the list of rivers and indices for this grid
    gri_fn = Ldir['grid'] / 'triv_info.csv'
    gri_df = pd.read_csv(gri_fn, index_col='rname')
    if Ldir['testing']:
        gri_df = gri_df.loc[['Birch Bay', 'Purdy Cr', 'Burley Cr', 'Perry Cr','McLane Cr'],:]

    # get list of overlapping rivers
    overlapping_trivs = gri_df[gri_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
    # consolidate overlapping rivers
    combined_names = trapsfun.combine_adjacent(overlapping_trivs)
    gri_df_no_ovrlp = pd.DataFrame(columns=gri_df.columns) 
    gri_df_no_ovrlp.index.name='rname'
    for trname in gri_df.index: # loop through original dataframe
        if trname in overlapping_trivs: # look for rivers that are in the list of duplicates
            name_index = np.where(overlapping_trivs == trname)[0][0] # get index in the list of duplicates
            if name_index%2 == 0: # even index means first occurence of duplicate
                newname = combined_names[int(name_index/2)] # combine names of duplicates
                gri_df_no_ovrlp.loc[newname] = gri_df.loc[trname] # add combined source to dataframe
            # Note: second duplicate will be dropped (so idir, isign, and uv will come from the first duplicate)
        else:
            gri_df_no_ovrlp.loc[trname] = gri_df.loc[trname] # if not a duplicate, then just copy over original info

    NTRIV = len(gri_df_no_ovrlp)
    # NTRIV = len(gri_df)

    # get the flow, temperature, and nutrient data for these days
    qtbio_triv_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type)

    # Add time coordinate
    triv_ds['river_time'] = (('river_time',), ot_vec)
    triv_ds['river_time'].attrs['units'] = Lfun.roms_time_units
    triv_ds['river_time'].attrs['long_name'] = 'river time'

    # Add river coordinate
    triv_ds['river'] = (('river',), np.arange(NRIV+1,NRIV+NTRIV+1))
    triv_ds['river'].attrs['long_name'] = 'tiny river runoff identification number'

    # Add Vshape
    vn = 'river_Vshape'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = ('s_rho', 'river')
    # For Vtransform = 2, even spacing is a good approximation, and
    # we implement this by using 1/N as the fraction in each vertical cell.
    Vshape = (1/N) * np.ones((N, NTRIV))
    triv_ds[vn] = (dims, Vshape)
    triv_ds[vn].attrs['long_name'] = vinfo['long_name']

    # Add position and direction
    for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        if vn == 'river_direction':
            triv_ds[vn] = (('river',), gri_df_no_ovrlp.idir.to_numpy().astype(int))
        elif vn == 'river_Xposition':
            X_vec = np.nan * np.ones(NTRIV)
            ii = 0
            for rn in gri_df_no_ovrlp.index:
                if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                    X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py'] + 1
                elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                    X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py']
                ii += 1
            triv_ds[vn] = (('river',), X_vec)
        elif vn == 'river_Eposition':
            E_vec = np.nan * np.ones(NTRIV)
            ii = 0
            for rn in gri_df_no_ovrlp.index:
                if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                    E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py']
                elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                    E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py'] + 1
                ii += 1
            triv_ds[vn] = (('river',), E_vec)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']

    # Add transport
    vn = 'river_transport'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('river',)
    Q_mat = np.zeros((NT, NTRIV))
    rr = 0
    for rn in gri_df_no_ovrlp.index:
        # sum flowrates together if duplicate river
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
        rr += 1
    triv_ds[vn] = (dims, Q_mat)
    triv_ds[vn].attrs['long_name'] = vinfo['long_name']
    triv_ds[vn].attrs['units'] = vinfo['units']

    # Add salinity and temperature
    for vn in ['river_salt', 'river_temp']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('s_rho', 'river')
        if vn == 'river_salt':
            TS_mat = np.zeros((NT, N, NTRIV))
        elif vn == 'river_temp':
            TS_mat = np.nan * np.zeros((NT, N, NTRIV))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                if '+' in rn:
                    # split into individual rivers
                    [riv1,riv2] = rn.split('+')
                    # get individual river dataframe
                    qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                    qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                    # calculate weighted average
                    temps = trapsfun.weighted_average('temp',qtbio_triv_df_1, qtbio_triv_df_2)
                else:
                    qtbio_triv_df = qtbio_triv_df_dict[rn]
                    temps = qtbio_triv_df['temp'].values
                for nn in range(N):
                    TS_mat[:, nn, rr] = temps
                rr += 1
        if np.isnan(TS_mat).any():
            print('Error from traps: nans in tiny river river_temp!')
            sys.exit()
        triv_ds[vn] = (dims, TS_mat)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']
        triv_ds[vn].attrs['units'] = vinfo['units']

    # Add biology that have existing climatology
    for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
        vn = 'river_' + var
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('s_rho', 'river')
        B_mat = np.nan * np.zeros((NT, N, NTRIV))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            if '+' in rn:
                # split into individual rivers
                [riv1,riv2] = rn.split('+')
                # get individual river dataframe
                qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                # calculate weighted average
                bvals = trapsfun.weighted_average(var,qtbio_triv_df_1, qtbio_triv_df_2)
            else:
                qtbio_triv_df = qtbio_triv_df_dict[rn]
                bvals = qtbio_triv_df[var].values
            for nn in range(N):
                B_mat[:, nn, rr] = bvals
            rr += 1
        if np.isnan(TS_mat).any():
            print('Error from traps: nans in tiny river bio!')
            sys.exit()
        triv_ds[vn] = (dims, B_mat)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']
        triv_ds[vn].attrs['units'] = vinfo['units']

    # Add remaining biology (see the lineup near the end of fennel_var.h)
    # I'm pretty sure this is simply filling everything with zeros
    bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
    for bvn in bvn_list:
        vn = 'river_' + bvn
        vinfo = zrfun.get_varinfo(vn)
        dims = (vinfo['time'],) + ('s_rho', 'river')
        B_mat = np.nan * np.zeros((NT, N, NTRIV))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            # qtbio_triv_df = qtbio_triv_df_dict[rn]
            for nn in range(N):
                B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
            rr += 1
        if np.isnan(B_mat).any():
            print('Error from traps: nans in B_mat for tiny river ' + vn)
            sys.exit()
        triv_ds[vn] = (dims, B_mat)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']
        triv_ds[vn].attrs['units'] = vinfo['units']

    # Rename rivers that share name with WWTP. This code appends ' R' at the end of the river name
    duplicates = ['Port Angeles', 'Port Townsend', 'Birch Bay', 'Port Gamble', 'Gig Harbor']
    print(gri_df_no_ovrlp.index)
    gri_df_no_ovrlp.index = np.where(gri_df_no_ovrlp.index.isin(duplicates), gri_df_no_ovrlp.index + ' R', gri_df_no_ovrlp.index)
    print(gri_df_no_ovrlp.index) 

    # Add river names
    triv_ds['river_name'] = (('river',), list(gri_df_no_ovrlp.index))
    triv_ds['river_name'].attrs['long_name'] = 'tiny river name'

###########################################################################################
# FORCING FOR MARINE POINT SOURCES

# Start Dataset
wwtp_ds = xr.Dataset()

NWWTP = len(gri_df)

if enable_pointsources == True:

    # Run placement algorithm to put marine point sources on LiveOcean grid
    trapsfun.traps_placement('wwtp')
    traps_type = 'wwtp'

    # define directory for tiny river climatology
    wwtp_dir = Ldir['LOo'] / 'pre' / 'traps' / 'point_sources'

    # climatological data files
    year0 = 1999
    year1 = 2017
    # climatological data
    Ldir['Cflow_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['Ctemp_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CDO_wwtp_fn']   = wwtp_dir / 'Data_historical' / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CNH4_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CNO3_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CTalk_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p')
    Ldir['CTIC_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p')

    # get the list of point sources and indices for this grid
    gri_fn = Ldir['grid'] / 'wwtp_info.csv'
    gri_df = pd.read_csv(gri_fn, index_col='rname')
    if Ldir['testing']:
        gri_df = gri_df.loc[['West Point', 'Birch Bay', 'Tacoma Central', 'US Oil & Refining'],:]
    gri_df = gri_df.drop('Birch Bay') # Remove the Birch Bay treatment plant
    
    # get list of overlapping point sources
    overlapping_wwtps = gri_df[gri_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
    # consolidate overlapping point sources
    combined_names = trapsfun.combine_adjacent(overlapping_wwtps)
    gri_df_no_ovrlp = pd.DataFrame(columns=gri_df.columns) 
    gri_df_no_ovrlp.index.name='rname'
    for psname in gri_df.index: # loop through original dataframe
        if psname in overlapping_wwtps: # look for point sources that are in the list of duplicates
            name_index = np.where(overlapping_wwtps == psname)[0][0] # get index in the list of duplicates
            if name_index%2 == 0: # even index means first occurence of duplicate
                newname = combined_names[int(name_index/2)] # combine names of duplicates
                gri_df_no_ovrlp.loc[newname] = gri_df.loc[psname] # add combined source to dataframe
            # Note: second duplicate will be dropped
        else:
            gri_df_no_ovrlp.loc[psname] = gri_df.loc[psname] # if not a duplicate, then just copy over original info

    NWWTP = len(gri_df_no_ovrlp)
    # NWWTP = len(gri_df)

    # get the flow, temperature, and nutrient data for these days
    qtbio_wwtp_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type)

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

    # Add Vshape
    vn = 'river_Vshape'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = ('s_rho', 'river')
    # All discharge coming from the bottom layer
    Vshape = np.zeros((N, NWWTP))
    Vshape[0,:] = 1
    wwtp_ds[vn] = (dims, Vshape)
    wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

    # Add position and direction
    for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        if vn == 'river_direction':
            wwtp_ds[vn] = (('river',), gri_df_no_ovrlp.idir.to_numpy().astype(int))
        elif vn == 'river_Xposition':
            X_vec = np.nan * np.ones(NWWTP)
            ii = 0
            for rn in gri_df_no_ovrlp.index:
                if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                    X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py'] + 1
                elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                    X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py']
                ii += 1
            wwtp_ds[vn] = (('river',), X_vec)
        elif vn == 'river_Eposition':
            E_vec = np.nan * np.ones(NWWTP)
            ii = 0
            for rn in gri_df_no_ovrlp.index:
                if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                    E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py']
                elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                    E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py'] + 1
                ii += 1
            wwtp_ds[vn] = (('river',), E_vec)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']

    # Add transport
    vn = 'river_transport'
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('river',)
    Q_mat = np.zeros((NT, NWWTP))
    rr = 0
    for rn in gri_df_no_ovrlp.index:
        # sum flowrates together if duplicate point sources
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
        Q_mat[:,rr] = flow * gri_df_no_ovrlp.loc[rn, 'isign']
        rr += 1
    wwtp_ds[vn] = (dims, Q_mat)
    wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
    wwtp_ds[vn].attrs['units'] = vinfo['units']

    # Add salinity and temperature
    for vn in ['river_salt', 'river_temp']:
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('s_rho', 'river')
        if vn == 'river_salt':
            TS_mat = np.zeros((NT, N, NWWTP))
        elif vn == 'river_temp':
            TS_mat = np.nan * np.zeros((NT, N, NWWTP))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                if '+' in rn:
                    # split into individual point sources
                    [wwtp1,wwtp2] = rn.split('+')
                    # get individual point source dataframe
                    qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                    qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                    # calculate weighted average
                    temps = trapsfun.weighted_average('temp',qtbio_wwtp_df_1, qtbio_wwtp_df_2)
                else:
                    qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                    temps = qtbio_wwtp_df['temp'].values
                for nn in range(N):
                    TS_mat[:, nn, rr] = temps
                rr += 1
        if np.isnan(TS_mat).any():
            print('Error from traps: nans in point source river_temp!')
            sys.exit()
        wwtp_ds[vn] = (dims, TS_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

    # Add biology that have existing climatology
    for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
        vn = 'river_' + var
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('s_rho', 'river')
        B_mat = np.nan * np.zeros((NT, N, NWWTP))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            if '+' in rn:
                # split into individual point sources
                [wwtp1,wwtp2] = rn.split('+')
                # get individual point source dataframe
                qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                # calculate weighted average
                bvals = trapsfun.weighted_average(var,qtbio_wwtp_df_1, qtbio_wwtp_df_2)
            else:
                qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                bvals = qtbio_wwtp_df[var].values
            for nn in range(N):
                B_mat[:, nn, rr] = bvals
            rr += 1
        if np.isnan(TS_mat).any():
            print('Error from traps: nans in tiny river bio!')
            sys.exit()
        wwtp_ds[vn] = (dims, B_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

    # Add remaining biology (see the lineup near the end of fennel_var.h)
    # I'm pretty sure this is simply filling everything with zeros
    bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
    for bvn in bvn_list:
        vn = 'river_' + bvn
        vinfo = zrfun.get_varinfo(vn)
        dims = (vinfo['time'],) + ('s_rho', 'river')
        B_mat = np.nan * np.zeros((NT, N, NWWTP))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            # qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
            for nn in range(N):
                B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
            rr += 1
        if np.isnan(B_mat).any():
            print('Error from traps: nans in B_mat for tiny river ' + vn)
            sys.exit()
        wwtp_ds[vn] = (dims, B_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

###########################################################################################

# combine all forcing datasets
all_ds = xr.merge([LOriv_ds,triv_ds, wwtp_ds])
# print(all_ds.river_name)

# Save to NetCDF
all_ds.to_netcdf(out_fn)
all_ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
