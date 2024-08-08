"""
Functions from Marvin Lorenz for making multi-layer TEF calculations, for bulk_calc.
Modifed lightly by Parker MacCready.
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
    
def calc_bulk_values(s, thisQ_dict, vn_list, print_info=False, min_trans=1):
    """
    input
    s=salinity array (sedges)
    Integrated transport arrays vs. S for a given time
    min_trans = minimum transport to consider
    """    
    # use the find_extrema algorithm
    ind, minmax = find_extrema(thisQ_dict['q'], print_info=print_info)
    
    # compute dividing salinities
    smin=s[0]
    DS=s[1]-s[0]
    div_sal=[]
    i=0
    while i < len(ind): 
        div_sal.append(smin+DS*ind[i])
        i+=1
                
    #calculate transports etc.
    in_dict = dict()
    out_dict = dict()
    for vn in vn_list:
            in_dict[vn] = []
            out_dict[vn] = []
    index=[]
    i=0
    vn_list_short = [item for item in vn_list if item != 'q']
    while i < len(ind)-1:
        # compute the transports and sort to in and out
        Q_i=-(thisQ_dict['q'][ind[i+1]]-thisQ_dict['q'][ind[i]])
        if Q_i<0 and np.abs(Q_i)>min_trans:
            out_dict['q'].append(Q_i)
        elif Q_i > 0 and np.abs(Q_i)>min_trans:
            in_dict['q'].append(Q_i)
        else:
            index.append(i)
            # this allows clipping of cases with tiny transport
        for vn in vn_list_short:
            F_i =-(thisQ_dict[vn][ind[i+1]] - thisQ_dict[vn][ind[i]])
            # vn_i=np.abs(F_i)/np.abs(Q_i)
            vn_i=(F_i)/(Q_i) # do not use absolute value (some mombal terms can be negative) Erin Broatch 2024.08.08
            if Q_i<0 and np.abs(Q_i)>1:
                out_dict[vn].append(vn_i)
            elif Q_i > 0 and np.abs(Q_i)>1:
                in_dict[vn].append(vn_i)
        i+=1
    # remove clipped div_sal values
    div_sal = np.delete(div_sal, index)
        
    return (in_dict, out_dict, div_sal, ind, minmax)