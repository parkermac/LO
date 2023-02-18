"""
Code to calculate a TEF time series using Marvin Lorenz' multi-layer code.
Based on his code, and modified by PM.

PERFORMANCE: Took 34 minutes for a full year with cas6_c0, only salt.

To test on mac:
run bulk_calc -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True

And for a full year:

(this only has salt)
run bulk_calc -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd

from lo_tools import Lfun, zfun
import tef_fun_lorenz as tfl

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# import tef_fun

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.p')]
if Ldir['testing']:
    sect_list = ['jdf3.p']

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

    # load the data file
    TEF = pickle.load(open(in_dir / snp, 'rb'))
    
    # add the absolute value of the net transport (to make Qprism)
    TEF['qabs'] = np.abs(TEF['qnet'].copy())
    
    # vn_list is variables that are arrays [ot, sbins]
    vn_list = [item for item in TEF.keys() if item not in ['sbins', 'ot', 'qnet', 'qabs', 'fnet', 'ssh']]
    
    # vec_list is time series [ot]
    vec_list = ['qnet', 'qabs', 'fnet', 'ssh']
    sbins = TEF['sbins']
    ot = TEF['ot']
    
    # do a little massaging of ot
    dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
    dt = dti.to_pydatetime() # an array of datetimes
    ot = np.array([Lfun.datetime_to_modtime(item) for item in dt])

    # tidal averaging, subsample, and cut off nans
    pad = 36
    # this pad is more than is required for the nans from the godin filter (35),
    # but, when combined with the subsampling we end up with fields at Noon of
    # each day (excluding the first and last days of the record)
    TEF_lp = dict()
    for vn in vn_list:
        #TEF_lp[vn] = zfun.filt_godin_mat(TEF[vn])[pad:-pad+1:24, :]
        TEF_lp[vn] = zfun.lowpass(TEF[vn], f='godin')[pad:-pad+1:24, :]
    for vn in vec_list:
        #TEF_lp[vn] = zfun.filt_godin(TEF[vn])[pad:-pad+1:24]
        TEF_lp[vn] = zfun.lowpass(TEF[vn], f='godin')[pad:-pad+1:24]
    ot = ot[pad:-pad+1:24]
    
    # also make an array of datetimes to save as the ot variable
    otdt = np.array([Lfun.modtime_to_datetime(item) for item in ot])
    
    if Ldir['testing']:
        print(Lfun.modtime_to_datetime(ot[0]))

    # get sizes and make sedges (the edges of sbins)
    DS=sbins[1]-sbins[0]
    sedges = np.concatenate((sbins,np.array([sbins[-1]] + DS))) - DS/2
    NT = len(ot)
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
                
    MLO['ot'] = otdt
    for vn in vec_list:
        MLO[vn] = TEF_lp[vn].copy()
            
    pickle.dump(MLO, open(out_fn, 'wb'))
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()
    
    if Ldir['testing']:
        plt.show()
        
print('\nTotal elapsed time = %d seconds' % (time()-tt00))





