# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:11:16 2016

@author: PM5

Some extra river functions.

"""


import os
import netCDF4 as nc
import pandas as pd
import river_class
import numpy as np
from datetime import datetime
import Lfun # assume path is provided by calling function

ncformat = 'NETCDF3_64BIT_OFFSET' # NETCDF3_CLASSIC'

def get_tc_rn(df):
    # makes a new column in df called 'tc_rn' which is the name of a
    # river that has T climatology
    for rn in df.index:
        if rn in ['coquille']:
            tc_rn = 'umpqua'
        elif rn in ['alsea']:
            tc_rn = 'siuslaw'
        elif rn in ['wilson', 'naselle', 'willapa', 'chehalis', 'humptulips',
                    'quinault', 'queets', 'hoh', 'calawah', 'hoko', 'elwha',
                    'dungeness', 'gold', 'sarita', 'sanjuan']:
            tc_rn = 'nehalem'
        elif rn in ['dosewallips', 'duckabush', 'hamma', 'skokomish', 'deschutes',
                    'nisqually', 'puyallup', 'green', 'snohomish',
                    'stillaguamish']:
            tc_rn = 'cedar'
        elif rn in ['skagit', 'samish']:
            tc_rn = 'nooksack'
        elif rn in ['clowhom', 'squamish']:
            tc_rn = 'fraser'
        elif rn in ['oyster', 'tsolum', 'englishman']:
            tc_rn = 'nanaimo'
        else:
            tc_rn = rn
        df.loc[rn, 'tc_rn'] = tc_rn
    return df
    
def get_qt(df, dt_ind, yd_ind, Ldir, dt1, days):
    
    dt_ind = dt_ind.tz_localize(tz=None)
    
    #%% step through all rivers
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    qt_df_dict = dict()
    for rn in df.index:
        print(10*'#' + ' ' + rn + ' ' + 10*'#')
        # initialize a qt (flow vs. time) DataFrame for this river
        qt_df = pd.DataFrame(index=dt_ind,
                             columns=['clim','his','usgs','ec','nws','final',
                                      'temperature'])
        rs = df.loc[rn] # a Series with info for this river
        riv = river_class.River(rn, rs)
        # Get climatology (using squeeze=True returns a Series)
        clim = pd.read_csv(str(Ldir['data']) + '/rivers/Data_clim/' + rn + '.csv',
                        header=None, index_col=0, squeeze=True)
        qt_clim_yd = clim.loc[yd_ind] # clip just the needed values
        # Get T climatology (using squeeze=True returns a Series)
        tc_rn = df.loc[rn, 'tc_rn']
        T_clim = pd.read_csv(str(Ldir['data']) + '/rivers/Data_T_clim/' + tc_rn + '.csv',
                        header=None, index_col=0, squeeze=True)
        T_clim_yd = T_clim.loc[yd_ind] # clip just the needed values
        # start to populate the qt DataFrame
        qt_df['clim'] = pd.Series(index=dt_ind, data=qt_clim_yd.values)

        # Get historical record (a Series)
        try:
            his = pd.read_pickle(str(Ldir['data']) + '/rivers/Data_historical/'
                        + rn + '.p')
            if dt1 <= his.index[-1]:
                # fill with historical data if the timing is right
                qt_df['his'] = his.reindex(dt_ind)
                qt_df['final'] = qt_df['his']
                print(' filled from historical')
            else:
                # otherwise try (sequentially) to fill from
                # nws, or usgs, or ec
                if pd.notnull(rs.nws) and Ldir['run_type'] == 'forecast':
                    riv.get_nws_data()
                    if not riv.qt.empty:
                        # debugging
                        #print(riv.qt)
                        riv.qt = riv.qt.tz_localize(tz=None)
                        qt_df['nws'] = riv.qt.reindex(dt_ind)
                        qt_df['final'] = qt_df['nws']
                        print(' filled from nws forecast')
                elif pd.notnull(rs.usgs):
                    riv.get_usgs_data(days)
                    if not riv.qt.empty:
                        # debugging
                        #print(riv.qt)
                        riv.qt = riv.qt.tz_localize(tz=None)
                        qt_df['usgs'] = riv.qt.reindex(dt_ind)
                        qt_df['final'] = qt_df['usgs']
                        print(' filled from usgs')
                        
                elif pd.notnull(rs.ec):
                    riv.get_ec_data(days)
                    if not riv.qt.empty:
                        # New 2019.02.14 to catch instances when
                        # we ask for data that is not in the last 18 months
                        # and the EC default is to return the most recent data
                        dt0_actual = riv.qt.index[0]
                        dt0_requested = days[0]
                        if np.abs((dt0_actual - dt0_requested).days) >= 1:
                            print(' request was out of range for ec')
                        else:
                            # debugging
                            # print(riv.qt)
                            # print(dt_ind)
                            riv.qt = riv.qt.tz_localize(tz=None)
                            qt_df['ec'] = riv.qt.reindex(dt_ind)
                            qt_df['final'] = qt_df['ec']
                            print(' filled from ec')
        except FileNotFoundError:
            # needed for analytical cases
            pass
        # check results and fill with extrapolation (ffill) or climatology
        # if False: # introduce errors for testing
        #     qt_df.loc[-3:, 'final'] = np.nan
        if ( pd.isnull(qt_df['final'].values).any() and
                not pd.isnull(qt_df['final'].values).all() ):
            qt_df['final'] = qt_df['final'].ffill(axis=0)
            print(' extended by ffill')
        if pd.isnull(qt_df['final'].values).any():
            qt_df['final'] = qt_df['clim']
            print( 'WARNING: missing values: all filled with climatology')
        if (qt_df['final'].values < 0).any():
            qt_df['final'] = qt_df['clim']
            print( 'WARNING: negative values: all filled with climatology')
        if pd.isnull(qt_df['final'].values).any():
            print( '>>>>>>> flow has missing values!! <<<<<<<<<')
        # Temperature data
        qt_df['temperature'] = pd.Series(index=dt_ind, data=T_clim_yd.values)
        if pd.isnull(qt_df['temperature'].values).any():
            print( '>>>>>>> temp has missing values!! <<<<<<<<<')
        # save in the dict
        qt_df_dict[rn] = qt_df
        
    return qt_df_dict
    
def dt64_to_dt(dt64):
    # convert numpy datetime64 to datetime
    dt = datetime.utcfromtimestamp(dt64.astype('datetime64[ns]').tolist()/1e9)
    return dt
    
def write_to_nc(out_fn, S, df, qt_df_dict, dt_ind):
    
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(out_fn, 'w', format=ncformat)
    
    nriv = len(df)
    N = S['N']
    ndt = len(dt_ind)
    SL = 50 # max string length for river names

    foo.createDimension('river', nriv)
    foo.createDimension('s_rho', N)
    foo.createDimension('river_time', ndt)
    foo.createDimension('slen', SL)

    v_var = foo.createVariable('river', float, ('river'))
    v_var[:] = np.arange(1, nriv+1)
    v_var.long_name = 'river runoff identification number'

    v_var = foo.createVariable('river_name', 'c', ('slen', 'river'))
    rr = 0
    for rn in df.index:
        print('- ' + rn)
        cc = 0
        for ch in rn:
            v_var[cc, rr] = ch
            cc += 1
        rr += 1

    v_var = foo.createVariable('river_time', float, ('river_time'))
    count = 0
    for item in dt_ind.values:
        item_dt = dt64_to_dt(item)
        v_var[count] = Lfun.datetime_to_modtime(item_dt)
        count += 1
    v_var.long_name = 'river runoff time'
    v_var.units = "seconds since 1970-01-01 00:00:00"

    v_var = foo.createVariable('river_direction', float, ('river'))
    count = 0
    for rn in df.index:
        v_var[count] = df.loc[rn, 'idir']
        count += 1
    v_var.long_name = 'river runoff direction'
    v_var.flag_values = "0, 1"
    v_var.flag_meanings = "flow across u-face, flow across v-face"
    v_varLwSrc_True = "flag not used"

    v_var = foo.createVariable('river_Xposition', float, ('river'))
    count = 0
    for rn in df.index:
        if df.loc[rn, 'idir'] == 0:
            v_var[count] = df.loc[rn, 'col_py'] + 1
        elif df.loc[rn, 'idir'] == 1:
            v_var[count] = df.loc[rn, 'col_py']
        count += 1
    v_var.long_name = 'river XI-position'
    v_var.LuvSrc_True_meaning = "i point index of U or V face source/sink"
    v_var.LwSrc_True_meaning = "i point index of RHO center source/sink" ;

    v_var = foo.createVariable('river_Eposition', float, ('river'))
    count = 0
    for rn in df.index:
        if df.loc[rn, 'idir'] == 0:
            v_var[count] = df.loc[rn, 'row_py']
        if df.loc[rn, 'idir'] == 1:
            v_var[count] = df.loc[rn, 'row_py'] + 1
        count += 1
    v_var.long_name = 'river ETA-position'
    v_var.LuvSrc_True_meaning = "j point index of U or V face source/sink"
    v_var.LwSrc_True_meaning = "j point index of RHO center source/sink" ;

    v_var = foo.createVariable('river_transport', float, ('river_time', 'river'))
    count = 0
    for rn in df.index:
        qt_df = qt_df_dict[rn]
        flow = qt_df['final'].values
        v_var[:, count] = flow * df.loc[rn, 'isign']
        count += 1
    v_var.long_name = 'river runoff vertically integrated mass transport'
    v_var.positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell"
    v_var.negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell"
    v_var.time = "river_time"
    v_var.units = "meter3 second-1"

    v_var = foo.createVariable('river_temp', float, ('river_time', 's_rho', 'river'))
    count = 0
    for rn in df.index:
        qt_df = qt_df_dict[rn]
        for nn in range(N):
            v_var[:, nn, count] = qt_df['temperature'].values
        count += 1
    v_var.long_name = 'river runoff potential temperature'
    v_var.time = "river_time"
    v_var.units = "Celsius"

    v_var = foo.createVariable('river_salt', float, ('river_time', 's_rho', 'river'))
    count = 0
    for rn in df.index:
        for nn in range(N):
            v_var[:, nn, count] = np.zeros(ndt)
        count += 1
    v_var.long_name = 'river runoff salinity'
    v_var.time = "river_time"
    v_var.units = "psu"
    
    # v_var = foo.createVariable('river_dye_01', float, ('river_time', 's_rho', 'river'))
    # count = 0
    # for rn in df.index:
    #     for nn in range(N):
    #         v_var[:, nn, count] = np.zeros(ndt)
    #     count += 1
    # v_var.long_name = 'river dye'
    # v_var.time = "river_time"
    # v_var.units = "kg m-3"

    v_var = foo.createVariable('river_Vshape', float, ('s_rho', 'river'))
    count = 0
    for rn in df.index:
        if False:
            # copied from old matlab code, and simplified
            #   %linear decay from surface value to 0, fixed by sng 7/2011
            v_var[:, count] = np.linspace(0,2/N,N)
            count += 1
        else:
            # new version 2019.01.20 PM, because the one above put
            # too much transport in the thin upper layers, which
            # caused VERY high velocities - and may have been the source
            # of blowups.
            csw = S['Cs_w']
            dcsw = np.diff(csw)
            v_var[:, count] = dcsw
            count += 1
            # should end up with velocity constant over depth
    v_var.long_name = 'river runoff mass transport vertical profile'
    v_var.requires = "must sum to 1 over s_rho"

    foo.close()
    
def add_bio(out_fn, df, yd_ind):
    vn_dict =  {'NO3':'MicroMolar',
               'phytoplankton':'MicroMolar N',
               'zooplankton':'MicroMolar N',
               'detritus':'MicroMolar N',
               'Ldetritus':'MicroMolar N',
               'CaCO3':'MicroMolar C',
               'oxygen':'MicroMolar O',
               'alkalinity':'MicroMolar',
               'TIC':'MicroMolar C'}
    foo = nc.Dataset(out_fn, 'a', format=ncformat)
    Vs = foo['river_Vshape'][:]
    N, nriv = Vs.shape
    for vn in vn_dict.keys():      
        v_var = foo.createVariable('river_'+vn, float, ('river_time', 's_rho', 'river'))
        count = 0
        for rn in df.index:
            vv = get_bio_vec(vn, rn, yd_ind)
            for nn in range(N):
                v_var[:, nn, count] = vv
            count += 1
        v_var.long_name = 'river runoff ' + vn
        v_var.time = "river_time"
        v_var.units = vn_dict[vn]
        
    foo.close()
    
def get_bio_vec(vn, rn, yd_ind):
    ndt = len(yd_ind)
    yd = yd_ind.values
    ovec = np.ones(ndt)
    if vn == 'NO3':
        if rn == 'fraser':
            vv = 2 + (13/2) + (13/2)*np.cos(2*np.pi*((yd-30)/366))
        elif rn == 'columbia':
            vv = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
        else:
            vv = 5 * ovec
    elif vn == 'oxygen':
        vv = 350 * ovec
    elif vn in ['alkalinity', 'TIC']:
        if rn in ['columbia', 'deschutes', 'duwamish']:
            vv = 1000 + ovec
        else:
            vv = 300 * ovec
    else:
        # all others fill with zeros
        vv = 0 * ovec
    return vv    
    
