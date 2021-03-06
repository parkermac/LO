"""
This is code for doing cast extractions.

Test on mac in ipython:
run extract_casts.py -gtx cas6_v3_lo8b -ro 2 -test True
run extract_casts.py -gtx cas6_v3_lo8b -ro 2 -test True -cruises newport_line

Run for real on perigee (after transferring the required sta_df.p)
python extract_casts.py -gtx cas6_v3_lo8b -ro 2 -cruises newport_line > np.log &
"""

from lo_tools import Lfun, zfun, zrfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

import pandas as pd
import subprocess
import xarray as xr
import numpy as np
from datetime import datetime

def get_cast(out_fn, fn, lon, lat):
    
    # This function does the cast extraction and saves it to a NetCDF file.
    
    # Find indices nearest to the location
    if Ldir['testing']:
        dt = datetime(2019,7,4,2)
        fng = get_his_fn_from_dt(Ldir, dt)
    else:
        fng = fn
        
    G, S, T = zrfun.get_basic_info(fng)
    Lon = G['lon_rho'][0,:]
    Lat = G['lat_rho'][:,0]
    
    # error checking
    if (lon < Lon[0]) or (lon > Lon[-1]):
        print('ERROR: lon out of bounds ' + out_fn.name)
        return
    if (lat < Lat[0]) or (lat > Lat[-1]):
        print('ERROR: lat out of bounds ' + out_fn.name)
        return

    ix = zfun.find_nearest_ind(Lon, lon)
    iy = zfun.find_nearest_ind(Lat, lat)
    
    # error checking
    if G['mask_rho'][iy,ix] == False:
        print('ERROR: point on land mask ' + out_fn.name)
        return
    
    if Ldir['testing'] == False:
        
        # Run ncks to do the extraction, overwriting any existing file
        cmd_list = ['ncks', '-d', 'xi_rho,'+str(ix), '-d', 'eta_rho,'+str(iy),
            '-v', 'AKs,salt,temp,NO3,phytoplankton,zooplankton,detritus,Ldetritus,oxygen,alkalinity,TIC,h',
            '-O', str(fn), str(out_fn)]
        # Note: We get AKs so that the s_w dimension is retained
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # and check on the results
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n' + ' sdtout '.center(60,'-'))
            print(stdout.decode())
        if len(stderr) > 0:
            print('\n' + ' stderr '.center(60,'-'))
            print(stderr.decode())
            
        # Add z-coordinates to the file using xarray
        foo = xr.load_dataset(out_fn)
        foo = foo.squeeze()
        z_rho, z_w = zrfun.get_z(foo['h'].values, np.array([0.]), S)
        foo['z_rho'] = (('s_rho'), z_rho)
        foo['z_w'] = (('s_w'), z_w)
        foo.s_rho.attrs['long_name'] = 'vertical position on s_rho grid, positive up, zero at surface'
        foo.s_rho.attrs['units'] = 'm'
        foo.s_w.attrs['long_name'] = 'vertical position on s_w grid, positive up, zero at surface'
        foo.s_w.attrs['units'] = 'm'
        foo.salt.attrs['units'] = 'g kg-1'
        foo.to_netcdf(out_fn)
        # and check on the results
        if Ldir['testing']:
            foo = xr.open_dataset(out_fn)
            for vn in foo.data_vars:
                print('%14s: %s' % (vn, str(foo[vn].shape)))
        foo.close()
    
    # plot for reality check
    if (Ldir['lo_env'] == 'pm_mac') and (Ldir['cruises'] == 'test_cruises'):
        a = xr.open_dataset(out_fn)
        import matplotlib.pyplot as plt
        plt.close('all')
        fig = plt.figure(figsize=(14,8))
        
        v_list = ['salt', 'temp', 'oxygen']
        NV = len(v_list)
        ii = 1
        for vn in v_list:
            ax = fig.add_subplot(1,NV,ii)
            ax.plot(a[vn].values, a['z_rho'].values, '-o')
            if ii == 1:
                ax.set_ylabel('Z [m]')
            ax.set_title('%s [%s]' % (vn, a[vn].units))
            ii += 1
        a.close()
        plt.show()
        
def get_his_fn_from_dt(Ldir, dt):
    # This creates the Path of a history file from its datetime
    date_string = dt.strftime(Ldir['ds_fmt'])
    his_num = ('0000' + str(dt.hour + 1))[-4:]
    fn = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + date_string) / ('ocean_his_' + his_num + '.nc')
    return fn

out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'cast' / Ldir['cruises']
Lfun.make_dir(out_dir, clean=True)

# Here is where we do the function call to actually make the cast extractions(s)
if Ldir['cruises'] == 'test_cruises':
    cruise = 'MyCruise'
    sn = 0
    out_fn = out_dir / (cruise + '_' + str(int(sn)) + '.nc')
    dt = datetime(2019,7,4,2)
    fn = get_his_fn_from_dt(Ldir, dt)
    lon = -123.228000
    lat = 48.240300
    get_cast(out_fn, fn, lon, lat)
    
elif Ldir['cruises'] == 'woac':
    sta_fn = Ldir['parent'] / 'ptools_output' / 'woac2' / 'sta_df.p'
    if sta_fn.is_file():
        sta_df = pd.read_pickle(sta_fn)
    cruises = sta_df['Cruise'].unique()
    for cruise in cruises:
        c_df = sta_df[sta_df['Cruise']==cruise]
        for sn in c_df.index:
            lon = c_df.loc[sn,'Longitude']
            lat = c_df.loc[sn,'Latitude']
            dt = c_df.loc[sn,'Datetime']
            out_fn = out_dir / (cruise + '_' + str(int(sn)) + '.nc')
            fn = get_his_fn_from_dt(Ldir, dt)
            if Ldir['testing']:
                print('\ncruise=%s sn=%d lon=%0.2f lat=%0.2f' % (cruise, sn, lon, lat))
                print(str(out_fn))
                print(fn)
            else:
                get_cast(out_fn, fn, lon, lat)
                
elif Ldir['cruises'] == 'wcoa':
    sta_fn = Ldir['parent'] / 'ptools_output' / 'wcoa' / 'sta_df.p'
    if sta_fn.is_file():
        sta_df = pd.read_pickle(sta_fn)
    cruises = sta_df['Cruise'].unique()
    for cruise in cruises:
        c_df = sta_df[sta_df['Cruise']==cruise]
        for sn in c_df.index:
            lon = c_df.loc[sn,'Longitude']
            lat = c_df.loc[sn,'Latitude']
            dt = datetime.strptime(c_df.loc[sn,'Datetime'],'%m/%d/%Y %H:%M:%S')
            out_fn = out_dir / (cruise + '_' + str(sn) + '.nc')
            fn = get_his_fn_from_dt(Ldir, dt)
            if Ldir['testing']:
                print('\ncruise=%s sn=%d lon=%0.2f lat=%0.2f' % (cruise, sn, lon, lat))
                print(str(out_fn))
                print(fn)
            else:
                get_cast(out_fn, fn, lon, lat)
                
elif Ldir['cruises'] == 'newport_line':
    sta_fn = Ldir['LOo'] / 'obs' / 'newport_line' / 'sta_df.p'
    if sta_fn.is_file():
        sta_df = pd.read_pickle(sta_fn)
    cruises = sta_df['Cruise'].unique()
    for cruise in cruises:
        c_df = sta_df[sta_df['Cruise']==cruise]
        for index, row in sta_df.iterrows():
            sn = row['Station']
            lon = row['Longitude']
            lat = row['Latitude']
            dt = row['Datetime']
            out_fn = out_dir / (cruise + '_' + str(sn) + '_' + dt.strftime(Lfun.ds_fmt) + '.nc')
            fn = get_his_fn_from_dt(Ldir, dt)
            get_cast(out_fn, fn, lon, lat)

