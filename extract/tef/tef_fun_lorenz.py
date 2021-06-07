# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:01:15 2018

@author: lorenz, edited by PM
"""

import numpy as np

def find_extrema(x, comp=5, print_info=False): # new, reduced, better described
    """
    input
    x = Q(S)
    comp = size of the window as an integer number
    """

    indices = []
    minmax = []
    
    i = 0
    while i < len(x): # np.shape(x)[0]:
        
        # set a to be i-1 (except for the first point, where a = 0)
        # assuming comp = 1
        if i-comp < 0:
            a = 0
        else:
            a=i-comp
            
        # set b to be i + 2 (unless it puts it out of bounds, then b = None)
        if i+comp >= len(x):
            b=None
        else:
            b=i+comp+1
            
        # so [a:b] is the three points centered on i (for comp = 1)
        if (x[i] == np.max(x[a:b])) and (np.max(x[a:b]) != np.min(x[a:b])):# and (x[i] != x[i-1]):
            indices.append(i)
            minmax.append('max')
        # this does not catch an initial increase...
        elif (x[i] == np.min(x[a:b])) and (np.max(x[a:b]) != np.min(x[a:b])):# and (x[i] != x[i-1]):
            indices.append(i)
            minmax.append('min')
        i+=1
        
    if print_info:
        print('* first step')
        print(' >> indices = %s' % (str(indices)))
        print(' >> minmax = %s' % (str(minmax)))
        
    # correct min min min or max max max parts,
    # especially in the beginning and end of the Q(S)
    ii=1
    while ii < len(indices):
        index=[]
        if minmax[ii] == minmax[ii-1]:
            if minmax[ii] == 'max': #note the index of the smaller maximum
                if x[indices[ii]]>=x[indices[ii-1]]:
                    index.append(ii-1)
                else:
                    index.append(ii)
            elif minmax[ii] == 'min': #note the index of the greater minimum
                if x[indices[ii]]<=x[indices[ii-1]]:
                    index.append(ii-1)
                else:
                    index.append(ii)
            minmax = np.asarray(minmax)
            indices = np.asarray(indices)
            indices = np.delete(indices, index)
            minmax = np.delete(minmax, index)
        else:
            ii+=1
    
    if print_info:
        print('* after correct min min min or max max max parts')
        print(' >> indices = %s' % (str(indices)))
        print(' >> minmax = %s' % (str(minmax)))
    
    # delete too small transports
    ii=0
    while ii < len(indices)-1: 
        index=[]
        
        # PM edit: define min_transport dynamically
        min_transport = (x.max() - x.min())/20
        if np.abs(x[indices[ii+1]]-x[indices[ii]]) <= min_transport:
            if ii == 0:
                # if smin is involved and the transport is too small,
                # smin has to change its min or max property
                index.append(ii+1)
                if minmax[ii] == 'min':
                    minmax[ii] = 'max'
                else:
                    minmax[ii] = 'min'
            elif ii+1==len(indices)-1:
                # if smax is involved and the transport is too small,
                # smin has to change its min or max property
                index.append(ii)
                if minmax[ii+1] == 'min':
                    minmax[ii+1] = 'max'
                else:
                    minmax[ii+1] = 'min'
            else: # else both involved div sals are kicked out
                if ii+2 < len(indices)-1:
                    # check and compare to i+2
                    if minmax[ii]=='min':
                        if x[indices[ii+2]]>x[indices[ii]]:
                            index.append(ii+2)
                            index.append(ii+1)
                        else:
                            index.append(ii)
                            index.append(ii+1)
                    elif minmax[ii]=='max':
                        if x[indices[ii+2]]<x[indices[ii]]:
                            index.append(ii+2)
                            index.append(ii+1)
                        else:
                            index.append(ii)
                            index.append(ii+1)
                else:
                    index.append(ii)
                    index.append(ii+1)
            # PM edit
            minmax = np.asarray(minmax)
            indices = np.asarray(indices)
            #
            indices = np.delete(indices, index)
            minmax = np.delete(minmax, index)
        else:
            ii+=1
            
    if print_info:
        print('* after delete too small transports')
        print(' >> indices = %s' % (str(indices)))
        print(' >> minmax = %s' % (str(minmax)))
    
    # so far the first and last minmax does not correspond
    # to smin and smax of the data, expecially smin due to numerical errors
    
    # correct smin index
    if len(x)>4: # an odd constraint - in practice this will always be true
        
        # this section seems to give rise to errors when the code above did not 
        # catch the initial min.
        ii=1
        while (np.abs(np.abs(x[ii])-np.abs(x[0])) < 1e-10) and (ii < len(x)-1):
            ii+=1
        indices[0]=ii-1
        #correct smax index
        if x[-1]==0: #for low salinity classes Q[-1] might not be zero as supposed.
            jj=-1
            while (x[jj] == 0) and (np.abs(jj) < len(x)-1):
                jj-=1
            indices[-1]=len(x)+jj+1
            
        # PM edit
        minmax = np.asarray(minmax)
        indices = np.asarray(indices)
        
    if print_info:
        print('* after correct smin index')
        print(' >> indices = %s' % (str(indices)))
        print(' >> minmax = %s' % (str(minmax)))
    
            
    return indices, minmax

def calc_bulk_values(s, Qv, Vv, Av, Qs, Qs2, print_info=False):
    """
    input
    s=salinity array
    Qv=Q(S)
    Qs=Q*s(S)
    Qs2=Q*s(S)*s(S)
    min_trans=minimum transport to consider
    """    
    # use the find_extrema algorithm
    ind, minmax = find_extrema(Qv, print_info=print_info)
    
    # compute dividing salinities
    smin=s[0]
    DS=s[1]-s[0]
    div_sal=[]
    i=0
    while i < len(ind): 
        div_sal.append(smin+DS*ind[i])
        i+=1
        
    #calculate transports etc.
    Q_in_m=[]
    Q_out_m=[]
    V_in_m=[]
    V_out_m=[]
    A_in_m=[]
    A_out_m=[]
    s_in_m=[]
    s_out_m=[]
    s2_in_m=[]
    s2_out_m=[]
    index=[]
    i=0
    while i < len(ind)-1:
        # compute the transports and sort to in and out
        Q_i=-(Qv[ind[i+1]]-Qv[ind[i]])
        V_i=Vv[ind[i]:ind[i+1]].mean()
        A_i=Av[ind[i]:ind[i+1]].sum()
        F_i =-( Qs[ind[i+1]] -Qs[ind[i]])
        F2_i=-(Qs2[ind[i+1]]-Qs2[ind[i]])
        #V_i=np.abs(Q_i)/np.abs(A_i) # doing this yields exactly the original TEF results
        s_i=np.abs(F_i)/np.abs(Q_i)
        s2_i=np.abs(F2_i)/np.abs(Q_i)
        if Q_i<0 and np.abs(Q_i)>1:
            Q_out_m.append(Q_i)
            V_out_m.append(V_i)
            A_out_m.append(A_i)
            s_out_m.append(s_i)
            s2_out_m.append(s2_i)
        elif Q_i > 0 and np.abs(Q_i)>1:
            Q_in_m.append(Q_i)
            V_in_m.append(V_i)
            A_in_m.append(A_i)
            s_in_m.append(s_i)
            s2_in_m.append(s2_i)
        else:
            index.append(i)
        i+=1
    div_sal = np.delete(div_sal, index)
        
    return Q_in_m, Q_out_m, V_in_m, V_out_m, A_in_m, A_out_m, s_in_m, s_out_m, s2_in_m, s2_out_m, div_sal, ind, minmax