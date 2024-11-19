"""
Code to experiment with downloading hycom fields and detiding them.
Apparently the new version has tides!!
"""

import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
from lo_tools import hycom_functions as hfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time

Ldir = Lfun.Lstart()

test_type = 'new_backfill'
verbose = True

# overide the spatial domain
hfun.aa = [-128, -127, 47, 48]

if test_type in ['original_forecast', 'new_forecast']:
    # Forecast version
    Ldir['date_string'] = '2024.11.19'
    this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
    nd_f = np.ceil(Ldir['forecast_days'])
    dt0 = this_dt - timedelta(days=2)
    dt1 = this_dt + timedelta(days=int(nd_f) + 2)
    dt_list_full = []
    dtff = dt0
    while dtff <= dt1:
        dt_list_full.append(dtff)
        dtff = dtff + timedelta(days=1)
elif test_type == 'new_backfill':
    # Several days of backfill
    Ldir['date_string'] = '2024.09.01'
    this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
    dt0 = this_dt
    dt1 = this_dt + timedelta(days=1)
    dt_list_full = []
    dtff = dt0
    while dtff <= dt1:
        dt_list_full.append(dtff)
        dtff = dtff + timedelta(days=1)

if test_type == 'original_forecast':
    # urls for extraction from new hycom (2024.09.14)
    url_ssh  = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ssh/FMRC_ESPC-D-V02_ssh_best.ncd'
    url_uvel_vvel = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_uv3z/FMRC_ESPC-D-V02_uv3z_best.ncd'
    url_temp_salt = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ts3z/FMRC_ESPC-D-V02_ts3z_best.ncd'
    # lists and dicts
    hkeys = ['ssh','vel','ts']
    url_list = [url_ssh, url_uvel_vvel, url_temp_salt]
    hycom_var_list = ["surf_el", "water_u,water_v", "water_temp,salinity"]
elif test_type == 'new_forecast':
    # New ones to try for forecast
    url_ssh  = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_ssh/FMRC_ESPC-D-V02_ssh_best.ncd'
    url_Sssh  = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_Sssh/FMRC_ESPC-D-V02_Sssh_best.ncd'
    url_uvel = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_u3z/FMRC_ESPC-D-V02_u3z_best.ncd'
    url_vvel = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_v3z/FMRC_ESPC-D-V02_v3z_best.ncd'
    url_temp = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_t3z/FMRC_ESPC-D-V02_t3z_best.ncd'
    url_salt = 'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_s3z/FMRC_ESPC-D-V02_s3z_best.ncd'
    # lists and dicts
    hkeys = ['Sssh','u', 'v','t', 's']
    url_list = [url_Sssh, url_uvel, url_vvel, url_temp, url_salt]
    hycom_var_list = ["steric_ssh", "water_u", "water_v", "water_temp", "salinity"]
elif test_type == 'new_backfill':
    # New ones to try for backfill
    url_ssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/ssh/2024'
    # url_Sssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/Sssh/2024'
    url_uvel  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/u3z/2024'
    url_vvel  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/v3z/2024'
    url_temp  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/t3z/2024'
    url_salt  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/s3z/2024'
    # lists and dicts
    hkeys = ['ssh','u']#, 'v','t', 's']
    url_list = [url_ssh, url_uvel]#, url_vvel, url_temp, url_salt]
    hycom_var_list = ['surf_el', "water_u"]#, "water_v", "water_temp", "salinity"]

url_dict = dict(zip(hkeys,url_list))
hycom_var_dict = dict(zip(hkeys,hycom_var_list))

h_out_dir = Ldir['LOo'] / 'forcing' / 'misc' / 'hycom_detide_test'
Lfun.make_dir(h_out_dir)

def messages(mess_str, stdout, stderr):
    # utility function to help with subprocess errors
    try:
        if len(stdout) > 0:
            print(mess_str)
            print(stdout.decode())
    except TypeError:
        pass
    try:
        if len(stderr) > 0:
            print(mess_str)
            print(stderr.decode())
    except TypeError:
        pass

def get_indices(h_out_dir, dt_list_full, verbose=False):
    # find and check the indices into hycom for the extraction

    # specify the sub region of hycom to extract
    aa = hfun.aa
    # convert to hycom format
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360

    ind_dicts = dict()
    for key in hkeys:
        got_indices = True

        if verbose:
            print('-- getting indices for: ' + key)
            sys.stdout.flush()
        url = url_dict[key]
        out_fn = h_out_dir / (key + '_tyx.nc')
        # get rid of the old version, if it exists
        out_fn.unlink(missing_ok=True)
        # extract coordinates
        cmd_list = ['ncks','-O','-v','time,lat,lon',url,str(out_fn)]
        #print(cmd_list)
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        messages('get_indices() messages:', stdout, stderr)

        if len(stderr) > 0:
            got_indices = False
            break

        # use the results
        ds = xr.open_dataset(out_fn)
        # find selected indices to use with ncks to extract fields
        t = ds.time.values
        tind = pd.DatetimeIndex(t)
        it_list = []
        for this_dt in dt_list_full:
            it = np.argwhere(tind==this_dt)[0][0]
            it_list.append(it)
        x = ds.lon.values
        y = ds.lat.values
        ix0 = zfun.find_nearest_ind(x,west)
        ix1 = zfun.find_nearest_ind(x,east)
        iy0 = zfun.find_nearest_ind(y,south)
        iy1 = zfun.find_nearest_ind(y,north)
        ds.close()
        ind_dict = {'it_list':it_list, 'ix0':ix0, 'ix1':ix1, 'iy0':iy0, 'iy1':iy1}
        ind_dicts[key] = ind_dict

    if got_indices:
        for key in hkeys:
            if ind_dicts[key] == ind_dict:
                pass
            else:
                print('Note from get_indices(): indices do not match for different variables.')
                """
                Looking at these the lon,lat indices were the same but the time indices were different for ssh.
                This is okay and we will deal with it in the extractions.
                """

    return ind_dicts, got_indices

def get_data_allday(idt, fld, out_fn, ind_dicts, verbose=False):
    """"
    This gets all the data in a day, which is hourly for ssh and 3-hourly for uvts.
    """
    # get rid of the old version, if it exists
    out_fn.unlink(missing_ok=True)
    print(' - getting hycom fields for ' + str(out_fn))
    got_fmrc = True
    
    tt0 = time()
    ind_dict = ind_dicts[fld]
    it0 = ind_dict['it_list'][idt]
    it1 = ind_dict['it_list'][idt+1]
    ix0 = ind_dict['ix0']
    ix1 = ind_dict['ix1']
    iy0 = ind_dict['iy0']
    iy1 = ind_dict['iy1']
    url = url_dict[fld]
    hvars = hycom_var_dict[fld]
    # extract data from HYCOM file
    AO = '-O'
    if fld in ['ssh','Sssh']:
        cmd_list = ['ncks',AO,'-v',hvars,
            '-d','time,'+str(it0)+','+str(it1),
            '-d','lat,'+str(iy0)+','+str(iy1),
            '-d','lon,'+str(ix0)+','+str(ix1),
            url,str(out_fn)]
    elif fld in ['u','v','t','s']:
        cmd_list = ['ncks',AO,'-v',hvars,
            '-d','time,'+str(it0)+','+str(it1),
            '-d','lat,'+str(iy0)+','+str(iy1),
            '-d','lon,'+str(ix0)+','+str(ix1),
            '-d','depth,'+str(0)+','+str(0),
            url,str(out_fn)]

    if False:
        print(cmd_list)
        sys.stdout.flush()
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    messages('ncks extract messages:', stdout, stderr)
            
    if len(stderr) > 0:
        got_fmrc = False

    print('Took %0.1f sec to get %s' % (time()-tt0, hvars))

    return got_fmrc

# Run the code
ind_dicts, got_indices = get_indices(h_out_dir, dt_list_full, verbose=verbose)

for idt in range(len(dt_list_full)-1):
    for fld in hkeys:
        got_fmrc = False # initialize each time
        data_out_fn =  h_out_dir / (fld + '_' + dt_list_full[idt].strftime(Lfun.ds_fmt)+ '.nc')
        if verbose:
            print('\n' + str(data_out_fn))
        sys.stdout.flush()
        got_fmrc = get_data_allday(idt, fld, data_out_fn, ind_dicts, verbose=verbose)


