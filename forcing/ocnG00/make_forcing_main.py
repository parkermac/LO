"""
This makes the ocn forcing files for ROMS, including the banas-fennel bio fields.

It uses GLORYS (glorys, hereafter) fields instead of HYCOM. Since the glorys fields
are daily averages is greatly simplifies the processing.

Testing:

run make_forcing_main.py -g cas7 -r backfill -d 2017.01.01 -f ocnG00 -test True


"""

from pathlib import Path
import sys
from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np
import pandas as pd

from lo_tools import Lfun, zfun, zrfun, Ofun2_nc
from lo_tools import glorys_functions as gfun
import Ofun
import Ofun_bio

# defaults
planB = False
add_CTD = False

# defaults related to testing
verbose = False
testing_planB = False

if Ldir['testing']:
    verbose = True
    from importlib import reload
    reload(gfun)
else:
    pass
    
# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# *** automate when to set add_CTD to True ***
if this_dt == datetime(2012,10,7):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True

# Set time limits for subsequent processing, or set planB to True.
if testing_planB == False:
    if Ldir['run_type'] == 'forecast':
        nd_f = np.ceil(Ldir['forecast_days'])
        dt0 = this_dt
        dt1 = this_dt + timedelta(days=int(nd_f))
        dt_list = pd.date_range(dt0,dt1,freq='D')
    elif Ldir['run_type'] == 'backfill':
        dt0 = this_dt
        dt1 = this_dt + timedelta(days=1)
    if verbose:
        print('dt0 = ' + str(dt0))
        print('dt1 = ' + str(dt1))
    dt_list = pd.date_range(dt0,dt1,freq='D')

else:
    # do this if testing_planB == True
    planB = True

# Get the required glorys extractions, or make sure we have them.
if planB == False:
    try:
        # set region to use
        if Ldir['gridname'] in ['cas7']:
            region = 'region7'
        else:
            print('unsupported region')
            sys.exit()
        # get data
        if Ldir['run_type'] == 'forecast':
            # Need to get fields from glorys site
            print('not supported yet')
            sys.exit()
        elif Ldir['run_type'] == 'backfill':
            # See if we already have the needed files
            # otherwise go to the glorys site and get them
            g_indir0 = Ldir['data'] / 'glorys' / region
            fn0 = g_indir0 / ('hindcast_' + dt0.strftime(Lfun.ds_fmt) + '.nc')
            fn1 = g_indir0 / ('hindcast_' + dt1.strftime(Lfun.ds_fmt) + '.nc')
            fn_list = [fn0, fn1]
            got_data = False
            if fn0.is_file() and fn1.is_file():
                got_data = True
            if got_data == True:
                planB = False
            elif got_data == False:
                print('- error getting data')
                planB = True
    except Exception as e:
        print(e)
        print('- error while getting glorys data')
        planB = True

if planB == False:
    try:
        # process the glorys fields to ROMS format

        V = dict() # dict to hold fields interpolated to ROMS grid

        G, S = gfun.create_roms_grid_info(Ldir['gridname'])

        # variables to look at
        vn_list = ['zeta','salt','temp','u','v']
        # vn_list = ['salt','temp','u','v']
        # vn_list = ['salt', 'zeta']

        # dict to relate the roms variables to the glorys variables
        vn_dict = {'salt':'so', 'temp':'thetao', 'u':'uo', 'v':'vo', 'zeta':'zos'}

        # NOTE: need to fill in ubar and vbar

        # get sizes
        NZ = S['N']; NR = G['M']; NC = G['L']

        # Make the time vector.
        ot_vec = np.array([Lfun.datetime_to_modtime(item) for item in dt_list])
        NT = len(ot_vec)

        # Create a dict of fields for the state variables.
        V = dict()
        VV = dict() # an extra dict for temporary things
        V['zeta'] = np.nan * np.zeros((NT, NR, NC))
        V['ubar'] = np.nan * np.zeros((NT, NR, NC-1))
        V['vbar'] = np.nan * np.zeros((NT, NR-1, NC))
        V['salt'] = np.nan * np.zeros((NT, NZ, NR, NC))
        V['temp'] = np.nan * np.zeros((NT, NZ, NR, NC))
        VV['z_rho'] = np.nan * np.zeros((NT, NZ, NR, NC))
        V['u'] = np.nan * np.zeros((NT, NZ, NR, NC-1))
        VV['z_u'] = np.nan * np.zeros((NT, NZ, NR, NC-1))
        VV['zw_u'] = np.nan * np.zeros((NT, NZ+1, NR, NC-1))
        V['v'] = np.nan * np.zeros((NT, NZ, NR-1, NC))
        VV['z_v'] = np.nan * np.zeros((NT, NZ, NR-1, NC))
        VV['zw_v'] = np.nan * np.zeros((NT, NZ+1, NR-1, NC))

        tt = 0 # time index counter
        for fng in fn_list:
            for vn in vn_list:
                vng = vn_dict[vn]
                print('Getting %s from %s' % (vn, vng))
                # if Ldir['run_type'] == 'forecast':
                #     if vng in ['so','thetao','zos']:
                #         fng = indir0 / 'glorys' / 'Data' / ('forecast_'+vng+'.nc')
                #     elif vng in ['uo','vo']:
                #         fng = indir0 / 'glorys' / 'Data' / 'forecast_cur.nc'
                # elif Ldir['run_type'] == 'backfill':
                #     fng = Ldir['parent'] / 'LPM_output' / 'glorys_archive' / ('hindcast_'+dstr+'.nc')
                if vn in ['salt', 'temp', 'zeta']:
                    gtag = 'rho'
                elif vn == 'u':
                    gtag = 'u'
                elif vn == 'v':
                    gtag = 'v'
                if vn == 'zeta':
                    do_2d = True
                else:
                    do_2d = False
                if do_2d == True:
                    FLD = gfun.interpolate_glorys_to_roms_2d(fng, vn, vng, gtag, G,
                        verbose=verbose, testing=Ldir['testing'])
                else:
                    zr, zw = gfun.get_zr(G['h'], V['zeta'][tt,:].squeeze(), S, vn)
                    if vn == 'salt':
                        z_rho = zr.copy()
                    elif vn == 'u':
                        z_u = zr.copy()
                        zw_u = zw.copy()
                    if vn == 'v':
                        z_v = zr.copy()
                        zw_v = zw.copy()
                    FLD = gfun.interpolate_glorys_to_roms(fng, vn, vng, gtag, zr, G,
                        verbose=verbose, testing=Ldir['testing'])
                V[vn][tt,:] = FLD
            VV['z_rho'][tt,:] = z_rho
            VV['z_u'][tt,:] = z_u
            VV['z_v'][tt,:] = z_v
            VV['zw_u'][tt,:] = zw_u
            VV['zw_v'][tt,:] = zw_v
            # note that the ":" here represents all axes after 0
            # and that it retains the correct shape
            tt += 1

    except Exception as e:
        print(e)
        print('- error processing glorys fields')
        planB = True

# Add extra fields
if planB == False:
    try:
        # create ubar and vbar
        for tt in range(NT):
            dzu = np.diff(VV['zw_u'][tt,:].squeeze(), axis=0)
            ubar = np.sum(V['u'][tt,:].squeeze() * dzu, axis=0) / np.sum(dzu, axis=0)
            dzv = np.diff(VV['zw_v'][tt,:].squeeze(), axis=0)
            vbar = np.sum(V['v'][tt,:].squeeze() * dzv, axis=0) / np.sum(dzv, axis=0)
            V['ubar'][tt,:] = ubar
            V['vbar'][tt,:] = vbar

        if add_CTD:
            tt0 = time()
            z_rho = zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
            for vn in ['salt','temp']:
                fld = V[vn].copy()
                V[vn] = Ofun_bio.fill_polygons(fld, vn, G, z_rho, Ldir)
            print(' - add_CTD task for salt and temp: %0.2f sec' % (time()-tt0))
        # add bio variables if needed
        if Ldir['do_bio']:
            tt0 = time()
            bvn_list = ['NO3', 'NH4', 'chlorophyll', 'phytoplankton', 'zooplankton',
                    'LdetritusN', 'SdetritusN', 'LdetritusC', 'SdetritusC',
                    'TIC', 'alkalinity', 'oxygen']
            salt = V['salt'].copy()
            for bvn in bvn_list:
                V[bvn] = Ofun_bio.create_bio_var(salt, bvn)
            print('- Add bio variables: %0.2f sec' % (time()-tt0))
            if add_CTD:
                tt0 = time()
                bvn_list_short = ['NO3','TIC', 'alkalinity', 'oxygen']
                z_rho = zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
                for bvn in bvn_list_short:
                    fld = V[bvn].copy()
                    V[bvn] = Ofun_bio.fill_polygons(fld, bvn, G, z_rho, Ldir)
                    V[bvn][mr3==0] = np.nan
                print(' - add_CTD task for bio variables: %0.2f sec' % (time()-tt0))
    except Exception as e:
        print(e)
        print('- error while adding extra fields to V')
        planB = True

# Write climatology file making use of zrfun.get_varinfo().
if planB == False:
    try:
        tt0 = time()
        out_fn = out_dir / 'ocean_clm.nc'
        out_fn.unlink(missing_ok=True)
        ds = xr.Dataset()    
        for vn in V.keys():
            # tt00 = time()
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            tname = vinfo['time_name']
            dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
            ds[vn] = (dims, V[vn])
            ds[vn].attrs['units'] = vinfo['units']
            ds[vn].attrs['long_name'] = vinfo['long_name']
            # time coordinate
            ds[tname] = ((tname,), ot_vec)
            ds[tname].attrs['units'] = Lfun.roms_time_units
        # add other coordinates
        for gtag in ['rho','u','v']:
            ds.coords['lon_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lon_'+gtag])
            ds.coords['lat_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lat_'+gtag])
        # add depth and z_rho
        ds['h'] = (('eta_rho', 'xi_rho'), G['h'])
        ds['z_rho_time'] = (('z_rho_time',), ot_vec)
        ds['z_rho_time'].attrs['units'] = Lfun.roms_time_units
        ds['z_rho'] = (('z_rho_time', 's_rho', 'eta_rho', 'xi_rho'), VV['z_rho'])
        # and save to NetCDF
        Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
        ds.to_netcdf(out_fn, encoding=Enc_dict)
        ds.close()
        print('- Write clm file: %0.2f sec' % (time()-tt0))
        sys.stdout.flush()
    except Exception as e:
        print(e)
        print('- error while creating clm file')
        planB = True

if True:

    if planB == True:
        # planB means that we use the ocean_clm.nc file from the day before and change its last time value
        # to be a day later.
        print('**** Using planB ****')
        result_dict['note'] = 'planB'
        ds_today = Ldir['date_string']
        dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
        dt_yesterday = dt_today - timedelta(days=1)
        ds_yesterday = datetime.strftime(dt_yesterday, format=Lfun.ds_fmt)
        clm_yesterday = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + ds_yesterday) / Ldir['frc'] / 'ocean_clm.nc'
        clm_today = out_dir / 'ocean_clm.nc'
        try:
            # use open_dataset, update, and save to a new name
            ds = xr.open_dataset(clm_yesterday, decode_times=False)
            tname_list = [item for item in ds.coords if 'time' in item]
            for tname in tname_list:
                ot_vec = ds[tname].values
                ot_vec[-1] += 86400
                ds.update({tname: (('ocean_time',), ot_vec)})
                ds[tname].attrs['units'] = Lfun.roms_time_units
            ds.to_netcdf(clm_today)
            ds.close()
        except Exception as e:
            print(e)

    if Ldir['start_type'] == 'new':
        # Write initial condition file if needed
        tt0 = time()
        in_fn = out_dir / 'ocean_clm.nc'
        out_fn = out_dir / 'ocean_ini.nc'
        out_fn.unlink(missing_ok=True)
        Ofun_nc.make_ini_file(in_fn, out_fn)
        print('- Write ini file: %0.2f sec' % (time()-tt0))
        sys.stdout.flush()

    # Write boundary file
    tt0 = time()
    in_fn = out_dir / 'ocean_clm.nc'
    out_fn = out_dir / 'ocean_bry.nc'
    out_fn.unlink(missing_ok=True)
    Ofun_nc.make_bry_file(in_fn, out_fn)
    print('- Write bry file: %0.2f sec' % (time()-tt0))
    sys.stdout.flush()

    def print_info(fn):
        print('\n' + str(fn))
        ds = xr.open_dataset(fn)#, decode_times=False)
        print(ds)
        ds.close()

    # Check results
    if Ldir['start_type'] == 'new':
        nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
    else:
        nc_list = ['ocean_clm.nc', 'ocean_bry.nc']
        
    if Ldir['testing']:
        # open datasets to have a peek manually
        dsc = xr.open_dataset(out_dir / 'ocean_clm.nc', decode_times=False)
        if Ldir['start_type'] == 'new':
            dsi = xr.open_dataset(out_dir / 'ocean_ini.nc', decode_times=False)
        dsb = xr.open_dataset(out_dir / 'ocean_bry.nc', decode_times=False)
            
    result_dict['result'] = 'SUCCESS'
    for fn in nc_list:
        if (out_dir / fn).is_file():
            pass
        else:
            result_dict['result'] = 'FAIL'

    # *******************************************************

    result_dict['end_dt'] = datetime.now()
    ffun.finale(Ldir, result_dict)
