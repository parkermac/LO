"""
This is the main program for making the ATM forcing file.

Test on mac in ipython:

short test:
run make_forcing_main.py -g cas6 -t v3 -r backfill -s continuation -d 2019.07.04 -test True -f atm0

test a forecast:
run make_forcing_main.py -g cas6 -t v3 -r forecast -s continuation -d 2019.07.04 -f atm0

test a forecast that will go to planB:
run make_forcing_main.py -g cas6 -t v3 -r forecast -s continuation -d 2019.07.05 -f atm0

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import os
import time
import shutil
import netCDF4 as nc
import numpy as np
import seawater as sw
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

from lo_tools import Lfun, zfun, zrfun

import atm_fun as afun
from importlib import reload
reload(afun)

# Set where are files located, and other situational choices.
do_d3 = True
do_d4 = True
wrf_dir = Ldir['data'] / 'wrf' # the default
if Ldir['lo_env'] == 'pm_mac':
    Ldir['run_type'] == 'backfill'
elif 'boiler' in Ldir['lo_env']:
    wrf_dir = Path('/data1/darr/wrf_crons/wrfout')
elif 'apogee' in Ldir['lo_env']:
    # wrf_dir = Path('/boildat/darr/wrf_crons/wrfout')
    wrf_dir = Path('/pgdat1/parker/LO_data/wrf')
elif 'perigee' in Ldir['lo_env']:
    wrf_dir = Path('/boildat1/darr/wrf_crons/wrfout')
# else:
#     print('WRF file location not yet supported on this machine.')
#     sys.exit()
    
# Create list of hours
if Ldir['run_type'] == 'backfill':
    hr_vec = range(25) # will iterate 0 to 24
elif Ldir['run_type'] == 'forecast':
    hr_max = Ldir['forecast_days'] * 24
    hr_vec = range(hr_max + 1)
    
# Create lists of input files.  These will be the complete lists
# regardless of whether or not the files exist, composed of string versions
# of the full path to each file.
d_str = Ldir['date_string'].replace('.','')
in_dir = wrf_dir / (d_str + '00')
d2_list = []
d3_list = []
d4_list = []
for hr in hr_vec:
    hr_str = ('0' + str(hr))[-2:]
    d2_list.append(in_dir / ('wrfout.ocean_d2.' + d_str + '00.f' + hr_str + '.0000'))
    d3_list.append(in_dir / ('wrfout.ocean_d3.' + d_str + '00.f' + hr_str + '.0000'))
    d4_list.append(in_dir / ('wrfout.ocean_d4.' + d_str + '00.f' + hr_str + '.0000'))
    
# Create dict that relates a d2 filename to a time index (used when writing to NetCDF)
d2i_dict = {}
for i, v in enumerate(d2_list):
    d2i_dict[v] = i
        
# Check for existence of files. If any d2 are missing then exit.
planB = False
for fn in d2_list:
    if not fn.is_file():
        print('** Missing file: ' + str(fn))
        planB = True
        break
    
if planB == False:
    # For d3 and d4 just make sure we have the first one, so that we can get the grid
    for fn in [d3_list[0]]:
        if not fn.is_file():
            print('** Missing file: ' + str(fn))
            do_d3 = False
    for fn in [d4_list[0]]:
        if not fn.is_file():
            print('** Missing file: ' + str(fn))
            do_d4 = False

    # Create vector of time, in model format (seconds since whenever)
    dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
    mod_time_list = []
    for hr in hr_vec:
        dt = dt0 + timedelta(days=hr/24)
        mod_time = Lfun.datetime_to_modtime(dt)
        mod_time_list.append(mod_time)
    mod_time_vec = np.array(mod_time_list)

    # Get ROMS model grid that we will interpolate to
    gds = nc.Dataset(Ldir['grid'] / 'grid.nc')
    lon = gds['lon_rho'][:]
    lat = gds['lat_rho'][:]
    gds.close()

    # Get WRF grids and grid size information (e.g. dx2_km = 12.5)
    # Note: lat, lon are only in the first file of the day (hour zero)
    lon2, lat2, dx2_km = afun.get_wrf_grid(d2_list[0])
    if do_d3:
        try:
            lon3, lat3, dx3_km = afun.get_wrf_grid(d3_list[0])
        except:
            do_d3 = False
    if do_d4:
        try:
            # sometimes there are empty files
            lon4, lat4, dx4_km = afun.get_wrf_grid(d4_list[0])
        except:
            do_d4 = False
    
    # Limit varlist if testing
    if Ldir['testing']:
        outvar_list = ['Pair']
    else:
        outvar_list = afun.outvar_list
    
    # Initialize NetCDF output files, one for each variable
    NR, NC = lon.shape
    NT = len(mod_time_list)
    nc_out_dict = {}
    for vn in outvar_list:
        # name output file
        out_fn = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc'] / (vn + '.nc')
        # print(out_fn)
        nc_out_dict[vn] = out_fn
        out_fn.unlink(missing_ok=True) # get rid of any old version
        foo = nc.Dataset(out_fn, 'w')
        # create dimensions
        timename = afun.timename_dict[vn]
        foo.createDimension(timename, NT) # could use None
        foo.createDimension('eta_rho', NR)
        foo.createDimension('xi_rho', NC)
        # add time data
        vv = foo.createVariable(timename, float, (timename,))
        vv.units = Lfun.roms_time_units
        vv[:] = mod_time_vec
        # add variable definition
        vv = foo.createVariable(vn, float, (timename, 'eta_rho', 'xi_rho'))
        vv.long_name = afun.longname_dict[vn]
        vv.units = afun.units_dict[vn]
        foo.close()

    # Find index to trim Eastern part of wrf fields
    lon_max = lon[0,-1] # easternmost edge of ROMS grid
    imax2 = zfun.find_nearest_ind(lon2[0,:], lon_max + .5)
    lon2 = lon2[:,:imax2]
    lat2 = lat2[:, :imax2]
    if do_d3:
        imax3 = zfun.find_nearest_ind(lon3[0,:], lon_max + .5)
        lon3 = lon3[:,:imax3]
        lat3 = lat3[:, :imax3]
    if do_d4:
        imax4 = zfun.find_nearest_ind(lon4[0,:], lon_max + .5)
        lon4 = lon4[:,:imax4]
        lat4 = lat4[:, :imax4]

    # Prepare coordinate arrays for interpolation
    XY = np.array((lon.flatten(), lat.flatten())).T # shape is (NR*NC, 2)
    XY2 = np.array((lon2.flatten(), lat2.flatten())).T
    if do_d3:
        XY3 = np.array((lon3.flatten(), lat3.flatten())).T
    if do_d4:
        XY4 = np.array((lon4.flatten(), lat4.flatten())).T

    # get nearest neighbor trees to use with wrf grids to interpolate
    # values from the wrf grids onto the ROMS grid
    IM2 = cKDTree(XY2).query(XY); IM2 = IM2[1]
    if do_d3:
        IM3 = cKDTree(XY3).query(XY); IM3 = IM3[1]
    if do_d4:
        IM4 = cKDTree(XY4).query(XY); IM4 = IM4[1]
    
    # Find coordinate rotation matrices to translate wrf velocity from
    # wrf grid directions to ROMS standard E+, N+
    ca2, sa2 = afun.get_angle(lon2, lat2)
    if do_d3:
        ca3, sa3 = afun.get_angle(lon3, lat3)
    if do_d4:
        ca4, sa4 = afun.get_angle(lon4, lat4)

    # Get the d3 and d4 masks. We appear to be avoiding a strip on the N.  Why?
    # We use these to overwrite just the parts where we have higher resolution wrf data.
    if do_d3:
        plon3_poly = np.concatenate((lon3[0,4:],lon3[:-5,-1],lon3[-5,4::-1],lon3[:-5:-1,4]))
        plat3_poly = np.concatenate((lat3[0,4:],lat3[:-5,-1],lat3[-5,4::-1],lat3[:-5:-1,4]))
        M3 = afun.get_indices_in_polygon(plon3_poly, plat3_poly, lon, lat)
    if do_d4:
        plon4_poly = np.concatenate((lon4[0,4:],lon4[:-5,-1],lon4[-5,4::-1],lon4[:-5:-1,4]))
        plat4_poly = np.concatenate((lat4[0,4:],lat4[:-5,-1],lat4[-5,4::-1],lat4[:-5:-1,4]))
        M4 = afun.get_indices_in_polygon(plon4_poly, plat4_poly, lon, lat)
    
    # MAIN TASK: loop over all hours

    if Ldir['testing']:
        # 20 = about noon local time
        d2_list = d2_list[20:22]
        d3_list = d3_list[20:22]
        d4_list = d4_list[20:22]
    dall_list = zip(d2_list, d3_list, d4_list)
    # Check out help(zip) to see how this works.  It creates an interable
    # that returns tuples made sequentially from entries of the things you zipped.
    # Note that this always works because we made our lists synthetically without regard
    # for if the files existed.
    for fn2, fn3, fn4 in dall_list:
        print('Working on ' + str(fn2).split('/')[-1] + ' and etc.')
    
        # flags to allow processing more files
        do_this_d3 = True
        do_this_d4 = True
    
        # if we are missing a d3 or d4 file then we don't work on it
        if not fn3.is_file():
            print(' - missing ' + str(fn3))
            do_this_d3 = False
        if not fn4.is_file():
            print(' - missing ' + str(fn4))
            do_this_d4 = False
    
        try:
            # On 2022.06.07 one of the d2 files showed up with a time index of length 0
            # (it should have been 1).  This is designed to respond to that specific error.
            # I could wrap the whole job in a try-except but I like to keep it more specific.
            ov2_dict = afun.gather_and_process_fields(fn2, imax2, ca2, sa2, outvar_list)
        except Exception as e:
            print('Error in gather and process d2 fields - going to Plan B')
            print(e)
            planB = True
            break
            
        ovi2_dict = afun.interp_to_roms(ov2_dict, outvar_list, IM2, NR, NC)
    
        if do_this_d3:
            try:
                ov3_dict = afun.gather_and_process_fields(fn3, imax3, ca3, sa3, outvar_list)
                ovi3_dict = afun.interp_to_roms(ov3_dict, outvar_list, IM3, NR, NC)
            except:
                print(' - could not process ' + str(fn3))
                do_this_d3 = False
    
        if do_this_d4:
            try:
                ov4_dict = afun.gather_and_process_fields(fn4, imax4, ca4, sa4, outvar_list)
                ovi4_dict = afun.interp_to_roms(ov4_dict, outvar_list, IM4, NR, NC)
            except:
                print(' - could not process ' + str(fn4))
                do_this_d4 = False
    
        # combine the grids
        ovc_dict = dict()
        for ovn in outvar_list:
            v2 = ovi2_dict[ovn]
            v = v2.copy()
            if do_this_d3:
                v3 = ovi3_dict[ovn]
                v[M3] = v3[M3]
            if do_this_d4:
                v4 = ovi4_dict[ovn]
                v[M4] = v4[M4]
            if np.sum(np.isnan(v)) > 0:
                print('** WARNING Nans in combined output ' + ovn)
            ovc_dict[ovn] = v
    
        # save to NetCDF
        tt = d2i_dict[fn2]
        for vn in outvar_list:
            fn = nc_out_dict[vn]
            foo = nc.Dataset(fn, 'a')
            foo[vn][tt,:,:] = ovc_dict[vn]
            foo.close()

if planB == True:
    # We make this an in instead of elif because planB might have been set to True
    # in the if planB == True section above.  Not great coding but I'm not sure how to
    # handle it more cleanly.
    result_dict['note'] = 'planB'
    print('**** Using planB ****')
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, Lfun.ds_fmt)
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format=Lfun.ds_fmt)
    
    LOogf_f_today = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']
    LOogf_f_yesterday = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + ds_yesterday) / Ldir['frc']
    
    outvar_list = afun.outvar_list
    nc_out_dict = dict()
    for ovn in outvar_list:
        fn_yesterday = LOogf_f_yesterday / (ovn + '.nc')
        fn_today = LOogf_f_today / (ovn + '.nc')
        nc_out_dict[ovn] = fn_today
        shutil.copyfile(fn_yesterday, fn_today)
        ds = nc.Dataset(fn_today, 'a')
        # advance the time by one day
        ot = ds[afun.timename_dict[ovn]][:]
        ot += 86400
        ds[afun.timename_dict[ovn]][:] = ot
        # and copy data from the previous backfill/forecast
        v = ds[ovn][:]
        print('%s %s' % (ovn, v.shape))
        NT, NR, NC = v.shape
        vv = v.copy()
        if NT == 25:
            # if it is backfill we just repeat the previous day
            pass
        elif NT == 73:
            # if it is forecast we use the last two days of the
            # previous forecast, and then repeat the last of these
            # as a best guess of day 3
            vv[:49,:,:] = v[24:,:,:]
            vv[49:,:,:] = v[49:,:,:]
        ds[ovn][:] = vv
        ds.close()
        
    if Ldir['lo_env'] == 'pm_mac':
        ds0 = nc.Dataset(fn_yesterday, 'r')
        ds1 = nc.Dataset(fn_today, 'r')
        ovn = 'Vwind'
        v0 = ds0[ovn][:,10,10]
        v1 = ds1[ovn][:,10,10]
        t0 = ds0[afun.timename_dict[ovn]][:]
        t1 = ds1[afun.timename_dict[ovn]][:]
        T0 = (t0-t0[0])/86400
        T1 = (t1-t0[0])/86400
        import matplotlib.pyplot as plt
        plt.close('all')
        plt.plot(T0,v0,'*r', T1,v1,'-k')
        plt.show()
        ds0.close()
        ds1.close()
        
# -------------------------------------------------------

# test for success
result_dict['result'] = 'success'
for vn in outvar_list:
    fn = nc_out_dict[vn]
    if fn.is_file():
        pass
    else:
       result_dict['result'] = 'fail'
       
# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
