"""
Code to calculate a TEF time series using Marvin Lorenz' multi-layer code.
Based on his code, and modified by PM.

PERFORMANCE: Takes about 5-27 minutes per year for 39 cas6 sections.

To test on mac:
run bulk_calc -gtagex cas6_v3_lo8b -0 2019.07.04 -1 2019.07.06

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time

from lo_tools import Lfun, zfun
import tef_fun_lorenz as tfl
from importlib import reload
reload(tfl)

Ldir = Lfun.Lstart()

# ---------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtagex', type=str, default='')   # e.g. cas6_v3_lo8b
parser.add_argument('-0', '--ds0', type=str, default='')        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='')        # e.g. 2019.07.04
args = parser.parse_args()

in_dir00 = Ldir['LOo'] / 'extract'
if len(args.gtagex) == 0:
    gtagex = Lfun.choose_item(in_dir00)
else:
    gtagex = args.gtagex
in_dir0 = in_dir00 / gtagex / 'tef'
if (len(args.ds0)==0) or (len(args.ds1)==0):
    ext_name = Lfun.choose_item(in_dir0, tag='processed')
else:
    ext_name = 'processed_' + args.ds0 + '_' + args.ds1
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.p')]
    
out_dir = in_dir0 / ext_name.replace('processed', 'bulk')
Lfun.make_dir(out_dir, clean=True)

# ---------

tt00 = time()

testing = False
# setting testing = True runs a deep debugging step, in which you only process
# the first day, and look at the details of the multi-layer bulk calculation,
# both graphically and as screen output.

# debugging
#sect_list = ['ai1.p', 'ss3.p']

for snp in sect_list:
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
    
    if testing:
        print(Lfun.modtime_to_datetime(ot.data[0]))

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
    
    if testing:
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
                
    MLO['ot'] = ot
    for vn in vec_list:
        MLO[vn] = TEF_lp[vn].copy()
        
    # debugging
    # print(MLO['ssh'].shape)
    # print(MLO['ssh'])
    
    pickle.dump(MLO, open(out_fn, 'wb'))
    
    if testing:
        plt.show()
        
print('\nTotal elapsed time = %d seconds' % (time()-tt00))





