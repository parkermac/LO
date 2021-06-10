# -*- coding: utf-8 -*-
"""
Code to calculate a TEF time series using Marvin Lorenz' new multi-layer code.

Based on his code, and modified by PM.

Takes about 10 minutes per year for 39 cas6 sections.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import matplotlib.pyplot as plt
import numpy as np
import pickle

import Lfun
import zfun

import tef_fun_lorenz as tfl
from importlib import reload
reload(tfl)

Ldir = Lfun.Lstart()

in_dir00 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir00)
in_dir0 = in_dir00 / gtagex / 'tef'
ext_name = Lfun.choose_item(in_dir0, tag='processed')
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.p')]
    
out_dir = in_dir0 / ext_name.replace('processed', 'bulk')
Lfun.make_dir(out_dir, clean=True)

# =================================

testing = True

for snp in sect_list:
    print('Working on ' + snp)
    out_fn = out_dir / snp

    # load the data file
    TEF = pickle.load(open(in_dir / snp, 'rb'))
    
    vn_list = [item for item in TEF.keys() if item not in ['sbins', 'ot', 'qnet', 'fnet', 'ssh']]
    vec_list = ['qnet', 'fnet', 'ssh']
    sbins = TEF['sbins']
    ot = TEF['ot']

    # tidal averaging, subsample, and cut off nans
    pad = 36
    # this pad is more than is required for the nans from the godin filter (35),
    # but, when combined with the subsampling we end up with fields at Noon of
    # each day (excluding the first and last days of the record)
    TEF_lp = dict()
    for vn in vn_list:
        TEF_lp[vn] = zfun.filt_godin_mat(TEF[vn])[pad:-pad+1:24, :]
    for vn in vec_list:
        TEF_lp[vn] = zfun.filt_godin(TEF[vn])[pad:-pad+1:24]
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
        MLO[vn] = TEF_lp[vn]
    pickle.dump(MLO, open(out_fn, 'wb'))
    
    if testing:
        plt.show()



