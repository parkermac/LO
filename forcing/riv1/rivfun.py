"""
Some extra river functions.

"""

import sys
import netCDF4 as nc
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from lo_tools import Lfun
from lo_tools import river_functions as rivf

def get_tc_rn(ri_df):
    """
    Makes a new column in df called 'tc_rn' which is the name of a
    river that has T climatology.
    
    As of 2021.04.19 we have a number of new rivers with T climatology,
    so this should be updated.  However it will work as is.
    """
    for rn in ri_df.index:
        if rn in ['coquille']:
            tc_rn = 'umpqua'
        elif rn in ['wilson']:
            tc_rn = 'nehalem'
        elif rn in ['naselle', 'willapa', 'chehalis', 'humptulips',
                    'queets', 'hoh', 'calawah', 'hoko', 'elwha',
                    'dungeness']:
            tc_rn = 'quinault'
        elif rn in ['dosewallips', 'duckabush', 'hamma', 'skokomish', 'deschutes',
                    'nisqually', 'puyallup', 'green', 'snohomish',
                    'stillaguamish']:#, 'nf_skokomish', 'sf_skokomish']:
            tc_rn = 'cedar'
        elif rn in ['samish']:
            tc_rn = 'nooksack'
        elif rn in ['squamish']:
            tc_rn = 'clowhom'
        elif rn in ['oyster', 'tsolum', 'englishman', 'cowichan']:
            tc_rn = 'nanaimo'
        elif rn in ['gold']:
            tc_rn = 'sarita'
        else:
            tc_rn = rn
        ri_df.loc[rn, 'tc_rn'] = tc_rn
    return ri_df
    
def get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days):
    # load historical and climatological data
    Hflow_df = pd.read_pickle(Ldir['Hflow_fn'])
    Cflow_df = pd.read_pickle(Ldir['Cflow_fn'])
    Ctemp_df = pd.read_pickle(Ldir['Ctemp_fn'])
    
    # initialize output dict
    qt_df_dict = dict()
    for rn in gri_df.index:
        
        # initialize screen output
        sout = '-- ' + rn + ': '
        
        rs = ri_df.loc[rn].copy() # a series with info for this river
        rs['got_data'] = False
        # initialize a qt (flow and temperature vs. time) DataFrame for this river
        qt_df = pd.DataFrame(index=dt_ind, columns=['clim','his','usgs','ec','nws','final','temperature'])
        # fill with historical and climatological fields
        qt_df.loc[:, 'clim'] = Cflow_df.loc[yd_ind,rn].values
        qt_df.loc[:, 'temperature'] = Ctemp_df.loc[yd_ind,ri_df.loc[rn,'tc_rn']].values
        try:
            qt_df.loc[:, 'his'] = Hflow_df.loc[dt_ind,rn].copy()
        except KeyError:
            pass # recent times will not be in the historical record

        if pd.notnull(qt_df.loc[:, 'his']).all():
            qt_df['final'] = qt_df['his']
            sout += 'filled from historical '
            rs.got_data = True
            
        # otherwise try (sequentially) to fill from
        # nws, or usgs, or ec
        if pd.notnull(rs.nws) and (Ldir['run_type'] == 'forecast'):
            rs, qt = rivf.get_nws_data(rs)
            if rs.got_data:
                qt_df['nws'] = qt.reindex(dt_ind)
                qt_df['final'] = qt_df['nws']
                sout += 'filled from nws forecast '
                
        if (not rs.got_data) and pd.notnull(rs.usgs):
            if rn in ['skokomish', 'hamma']:
                rs, qt = rivf.get_usgs_data_custom(rs, days)
            else:
                rs, qt = rivf.get_usgs_data(rs, days)
            if rs.got_data:
                qt_df['usgs'] = qt
                qt_df['final'] = qt_df['usgs']
                sout += 'filled from usgs '
                
        if (not rs.got_data) and pd.notnull(rs.ec):
            rs, qt = rivf.get_ec_data(rs, days)
            if rs.got_data:
                qt_df['ec'] = qt
                qt_df['final'] = qt_df['ec']
                sout += 'filled from ec '
        
        # check results and fill with extrapolation (ffill) or climatology
        if ( pd.isnull(qt_df['final'].values).any() and
                not pd.isnull(qt_df['final'].values).all() ):
            qt_df['final'] = qt_df['final'].ffill(axis=0)
            sout += ':: extended by ffill'
        if pd.isnull(qt_df['final'].values).any():
            qt_df['final'] = qt_df['clim']
            sout += '!! still missing values after ffill - all filled with climatology'
        if (qt_df['final'].values < 0).any():
            qt_df['final'] = qt_df['clim']
            sout += '!! negative values - all filled with climatology'
        if pd.isnull(qt_df['final'].values).any():
            sout += '>>>>>>> flow has missing values!! <<<<<<<<<'
        # Temperature data
        if pd.isnull(qt_df['temperature'].values).any():
            sout += '>>>>>>> temp has missing values!! <<<<<<<<<'
        # save in the dict
        qt_df_dict[rn] = qt_df
        
        # screen output
        print(sout)
        sys.stdout.flush()
        
    return qt_df_dict
        
def write_to_nc(out_fn, S, gri_df, qt_df_dict, dt_ind):
        
    out_fn.unlink(missing_ok=True)# get rid of the old version, if it exists
    foo = nc.Dataset(out_fn, 'w')
    
    nriv = len(gri_df.index)
    N = S['N']
    ndt = len(dt_ind)

    foo.createDimension('river', nriv)
    foo.createDimension('s_rho', N)
    foo.createDimension('river_time', ndt)

    v_var = foo.createVariable('river', float, ('river',))
    v_var[:] = np.arange(1, nriv+1)
    v_var.long_name = 'river runoff identification number'
    
    # make a dict to relate river name to a number (not the same as the river number above)
    rlist = gri_df.index.to_list()
    rc = dict(zip(rlist, range(nriv)))

    v_var = foo.createVariable('river_name', str, ('river',))
    for rn in gri_df.index:
        v_var[rc[rn]] = rn

    v_var = foo.createVariable('river_time', float, ('river_time',))
    count = 0
    for item in dt_ind:
        v_var[count] = Lfun.datetime_to_modtime(item)
        count += 1
    v_var.long_name = 'river runoff time'
        
    v_var.units = Lfun.roms_time_units

    v_var = foo.createVariable('river_direction', float, ('river',))
    for rn in gri_df.index:
        v_var[rc[rn]] = gri_df.loc[rn, 'idir']
    v_var.long_name = 'river runoff direction'
    v_var.flag_values = "0, 1"
    v_var.flag_meanings = "flow across u-face, flow across v-face"
    v_var.LwSrc_True = "flag not used"

    v_var = foo.createVariable('river_Xposition', float, ('river',))
    for rn in gri_df.index:
        if gri_df.loc[rn, 'idir'] == 0:
            v_var[rc[rn]] = gri_df.loc[rn, 'col_py'] + 1
        elif gri_df.loc[rn, 'idir'] == 1:
            v_var[rc[rn]] = gri_df.loc[rn, 'col_py']
    v_var.long_name = 'river XI-position'
    v_var.LuvSrc_True_meaning = "i point index of U or V face source/sink"
    v_var.LwSrc_True_meaning = "i point index of RHO center source/sink" ;

    v_var = foo.createVariable('river_Eposition', float, ('river',))
    for rn in gri_df.index:
        if gri_df.loc[rn, 'idir'] == 0:
            v_var[rc[rn]] = gri_df.loc[rn, 'row_py']
        if gri_df.loc[rn, 'idir'] == 1:
            v_var[rc[rn]] = gri_df.loc[rn, 'row_py'] + 1
    v_var.long_name = 'river ETA-position'
    v_var.LuvSrc_True_meaning = "j point index of U or V face source/sink"
    v_var.LwSrc_True_meaning = "j point index of RHO center source/sink" ;

    v_var = foo.createVariable('river_transport', float, ('river_time', 'river'))
    for rn in gri_df.index:
        qt_df = qt_df_dict[rn]
        flow = qt_df['final'].values
        v_var[:, rc[rn]] = flow * gri_df.loc[rn, 'isign']
    v_var.long_name = 'river runoff vertically integrated mass transport'
    v_var.positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell"
    v_var.negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell"
    v_var.time = "river_time"
    v_var.units = "meter3 second-1"

    v_var = foo.createVariable('river_temp', float, ('river_time', 's_rho', 'river'))
    for rn in gri_df.index:
        qt_df = qt_df_dict[rn]
        for nn in range(N):
            v_var[:, nn, rc[rn]] = qt_df['temperature'].values
    v_var.long_name = 'river runoff potential temperature'
    v_var.time = "river_time"
    v_var.units = "Celsius"

    v_var = foo.createVariable('river_salt', float, ('river_time', 's_rho', 'river'))
    for rn in gri_df.index:
        for nn in range(N):
            v_var[:, nn, rc[rn]] = np.zeros(ndt)
    v_var.long_name = 'river runoff salinity'
    v_var.time = "river_time"
    v_var.units = "psu"
    
    # v_var = foo.createVariable('river_dye_01', float, ('river_time', 's_rho', 'river'))
    # for rn in gri_df.index:
    #     for nn in range(N):
    #         v_var[:, nn, rc[rn]] = np.zeros(ndt)
    # v_var.long_name = 'river dye'
    # v_var.time = "river_time"
    # v_var.units = "kg m-3"

    v_var = foo.createVariable('river_Vshape', float, ('s_rho', 'river'))
    for rn in gri_df.index:
        # for Vtransform = 2 even spacing is a better approximation
        v_var[:, rc[rn]] = 1/N
        # goal is to end up with velocity constant over depth
    v_var.long_name = 'river runoff mass transport vertical profile'
    v_var.requires = "must sum to 1 over s_rho"

    foo.close()
    
def add_bio(out_fn, gri_df, yd_ind):
    
    # make a dict to relate river name to a number
    rlist = gri_df.index.to_list()
    rc = dict(zip(rlist, range(len(rlist))))
    
    vn_dict =  {'NO3':'MicroMolar',
               'phytoplankton':'MicroMolar N',
               'zooplankton':'MicroMolar N',
               'detritus':'MicroMolar N',
               'Ldetritus':'MicroMolar N',
               'CaCO3':'MicroMolar C',
               'oxygen':'MicroMolar O',
               'alkalinity':'MicroMolar',
               'TIC':'MicroMolar C'}
    foo = nc.Dataset(out_fn, 'a')
    Vs = foo['river_Vshape'][:]
    N, nriv = Vs.shape
    for vn in vn_dict.keys():      
        v_var = foo.createVariable('river_'+vn, float, ('river_time', 's_rho', 'river'))
        for rn in gri_df.index:
            vv = get_bio_vec(vn, rn, yd_ind)
            for nn in range(N):
                v_var[:, nn, rc[rn]] = vv
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
            vv = 1000 * ovec
        else:
            vv = 300 * ovec
    else:
        vv = 0 * ovec # all others filled with zeros
    return vv    
    
