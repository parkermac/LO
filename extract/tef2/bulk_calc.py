"""
Code to calculate a TEF time series using Marvin Lorenz' multi-layer code.
Based on his code, and modified by PM.

PERFORMANCE: 9 sec for test.

To test on mac:
run bulk_calc.py -gtx cas7_trapsV00_meV00 -ctag c0 -0 2017.07.04 -1 2017.07.06

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
import tef_fun_lorenz as tfl

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    sect_list = ['jdf3.nc']

# ---------

tt00 = time()

# setting Ldir['testing'] = True runs a deep debugging step, in which you only process
# the first day, and look at the details of the multi-layer bulk calculation,
# both graphically and as screen output.

for snp in sect_list:
    tt0 = time()
    
    print('Working on ' + snp)
    sys.stdout.flush()
    out_fn = out_dir / snp

    # load the processed Dataset for this section
    
    ds = xr.open_dataset(in_dir / snp)
    
    # Create the absolute value of the net transport (to make Qprism)
    # but first remove the low-passed transport (like Qr)
    qnet_lp = zfun.lowpass(ds.qnet.values, f='godin',nanpad=False)
    qabs = np.abs(ds.qnet.values - qnet_lp)
    
    # Tidal averaging, subsample, and cut off nans
    pad = 36
    # this pad is more than is required for the nans from the godin filter (35),
    # but, when combined with the subsampling we end up with fields at Noon of
    # each day (excluding the first and last days of the record)
    TEF_lp = dict() # temporary storage
    vn_list = []
    vec_list = []
    for vn in ds.data_vars:
        if ('time' in ds[vn].coords) and ('sbins' in ds[vn].coords):
            TEF_lp[vn] = zfun.lowpass(ds[vn].values, f='godin')[pad:-pad+1:24, :]
            vn_list.append(vn)
        elif ('time' in ds[vn].coords) and ('sbins'  not in ds[vn].coords):
            TEF_lp[vn] = zfun.lowpass(ds[vn].values, f='godin')[pad:-pad+1:24]
            vec_list.append(vn)
    time_lp = ds.time.values[pad:-pad+1:24]
    sbins = ds.sbins.values
    TEF_lp['qabs'] = zfun.lowpass(qabs, f='godin')[pad:-pad+1:24]
    # Add the qprism time series.
    # Conceptually, qprism is the maximum possible exchange flow if all
    # the flood tide made Qin and all the ebb tide made Qout.
    # If you go through the trigonometry you find that qprism = 1/2 <qabs>.
    TEF_lp['qprism'] = TEF_lp['qabs'].copy()/2
    vec_list += ['qabs', 'qprism']
    ds.close()
    
    # get sizes and make sedges (the edges of sbins)
    DS=sbins[1]-sbins[0]
    sedges = np.concatenate((sbins,np.array([sbins[-1]] + DS))) - DS/2
    NT = len(time_lp)
    NS = len(sedges)

    # calculate all transports integrated over salinity, e.g. Q(s) = integral(q ds)
    omat = np.zeros((NT, NS))
    Q_dict = dict()
    for vn in vn_list:
        Q_dict[vn] = omat.copy()
        Q_dict[vn][:,:-1] = np.fliplr(np.cumsum(np.fliplr(TEF_lp[vn]), axis=1))

    # prepare arrays to hold multi-layer output
    nlay = 30
    nmat = np.nan * np.ones((NT, nlay))
    MLO = dict()
    for vn in vn_list:
        MLO[vn] = nmat.copy()
    
    if Ldir['testing']:
        plt.close('all')
        dd_list = [0]
        print_info = True
    else:
        dd_list = range(NT)
        print_info = False

    for dd in dd_list:
            
        thisQ_dict = dict()
        for vn in vn_list:
            thisQ_dict[vn] = Q_dict[vn][dd,:]
    
        if print_info == True:
            print('\n**** dd = %d ***' % (dd))
                
        out_tup = tfl.calc_bulk_values(sedges, thisQ_dict, vn_list, print_info=print_info)
        in_dict, out_dict, div_sal, ind, minmax = out_tup
        
        if print_info == True:
            print(' ind = %s' % (str(ind)))
            print(' minmax = %s' % (str(minmax)))
            print(' div_sal = %s' % (str(div_sal)))
            print(' Q_in_m = %s' % (str(in_dict['q'])))
            print(' s_in_m = %s' % (str(in_dict['salt'])))
            print(' Q_out_m = %s' % (str(out_dict['q'])))
            print(' s_out_m = %s' % (str(out_dict['salt'])))
        
            fig = plt.figure(figsize=(12,8))
        
            ax = fig.add_subplot(121)
            ax.plot(Q_dict['q'][dd,:], sedges,'.k')
            min_mask = minmax=='min'
            max_mask = minmax=='max'
            print(min_mask)
            print(max_mask)
            ax.plot(Q_dict['q'][dd,ind[min_mask]], sedges[ind[min_mask]],'*b')
            ax.plot(Q_dict['q'][dd,ind[max_mask]], sedges[ind[max_mask]],'*r')
            ax.grid(True)
            ax.set_title('Q(s) Time index = %d' % (dd))
            ax.set_ylim(-.1,36.1)
            ax.set_ylabel('Salinity')
        
            ax = fig.add_subplot(122)
            ax.plot(TEF_lp['q'][dd,:], sbins)
            ax.grid(True)
            ax.set_title('-dQ/ds')
        
        bulk_dict = dict()
        for vn in vn_list:
            bulk_dict[vn] = np.array(in_dict[vn] + out_dict[vn])
        ii = np.argsort(bulk_dict['salt'])
        if len(ii) > 0:
            for vn in vn_list:
                bulk_dict[vn] = bulk_dict[vn][ii]
                NL = len(ii)
                MLO[vn][dd, :NL] = bulk_dict[vn]
                
    for vn in vec_list:
        MLO[vn] = TEF_lp[vn].copy()
        
    # Pack results in a Dataset and then save to NetCDF
    ds = xr.Dataset(coords={'time': time_lp,'layer': np.arange(nlay)})
    for vn in vn_list:
        ds[vn] = (('time','layer'), MLO[vn])
    for vn in vec_list:
        ds[vn] = (('time'), MLO[vn])
    # save it to NetCDF
    ds.to_netcdf(out_dir / out_fn)
    if Ldir['testing']:
        pass
    else:
        ds.close()
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()
    
    if Ldir['testing']:
        plt.show()
        
print('\nTotal elapsed time = %d seconds' % (time()-tt00))





