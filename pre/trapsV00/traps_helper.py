"""
helper functions for traps
"""
import numpy as np
import pandas as pd
from lo_tools import zfun

#################################################################################
#                          Make climatology helpers                             #
#################################################################################

def ds_to_avgdf(source_name,ecology_data_ds):
    '''
    Converts a dataset of Ecology's data (read from .nc file)
    into a dataframe for a single source. 
    Also outputs dataframe of the average year,
    and dataframe of the standard deviation over a year.
    '''
    d = {'Date': ecology_data_ds.date.values,
         'Flow(m3/s)':  ecology_data_ds.flow[ecology_data_ds.name==source_name,:].values[0],
         'Temp(C)':     ecology_data_ds.temp[ecology_data_ds.name==source_name,:].values[0],
         'NO3(mmol/m3)':ecology_data_ds.NO3[ecology_data_ds.name==source_name,:].values[0],
         'NH4(mmol/m3)':ecology_data_ds.NH4[ecology_data_ds.name==source_name,:].values[0],
         'TIC(mmol/m3)':ecology_data_ds.TIC[ecology_data_ds.name==source_name,:].values[0],
         'Talk(meq/m3)':ecology_data_ds.Talk[ecology_data_ds.name==source_name,:].values[0],
         'DO(mmol/m3)': ecology_data_ds.DO[ecology_data_ds.name==source_name,:].values[0]}
    df = pd.DataFrame(data=d)
    # replace all zeros with nans, so zeros don't bias data
    df = df.replace(0, np.nan)
    # add day of year column
    df['day_of_year'] = df.apply(lambda row: row.Date.dayofyear, axis = 1)
    # add year column
    df['year'] = pd.DatetimeIndex(df['Date']).year

    # calculate averages
    # (compress 1999-2017 timeseries to single year, with an average for each day)
    avgs_df = df.groupby('day_of_year').mean().reset_index()
    # calculate standard deviation
    sds_df = df.groupby('day_of_year').std(ddof=0).reset_index()

    # replace all nans with zeros, so I'm no longer injecting nans
    avgs_df = avgs_df.replace(np.nan,0)
    sds_df = sds_df.replace(np.nan,0)

    return df, avgs_df, sds_df

#################################################################################
#                          TRAPS placement helpers                              #
#################################################################################

def in_domain(x, y, X, Y):
    '''
    Note: code borrowed from pgrid/carve_rivers.py
    Utility function to make sure that a point (x, y) is
    in a domain specified by vectors X and Y.
    We actually require the point to be 'pad' in from the edge.

    Inputs:
        x:  point of interest lon coordinates
        y:  point of interest lat coordinates
        X:  domain lon coordinates
        Y:  domain lat coordinates
    Outputs:
        boolean (T if in domain, F if not in domain)
    '''
    pad = 1
    if x>=X[0+pad] and x<=X[-1-pad] and y>=Y[0+pad] and y<=Y[-1-pad]:
        return True
    else:
        return False

def cell_in_domain(ival,jval,II,JJ):
    '''
    Utility function to make sure that a grid cell with index location (ival,jval)
    in the domain specified by domain sizes II and JJ.

    Inputs:
        ival:   grid cell of interest i-index
        jval:   grid cell of interest j-index
        II:     domain max i index
        JJ:     domain max j index
    Outputs:
        boolean (T if in domain, F if not in domain)
    '''
    if ival>=0 and ival<=II and jval>=0 and jval<=JJ:
        return True
    else:
        return False

def get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho):
    '''
    Function that checks whether a grid cell is in the model domain and if it is on water.
    If both, then records the index of the grid cell, and the distance from
    the grid cell to the original WWTP lat/lon coordinates.

    Inputs:
        I_ind:      i-index of WWTP
        J_ind:      j-index of WWTP
        X:          domain lon coordinates
        Y:          domain lat coordinates
        x:          WWTP lon coordinates
        y:          WWTP lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        ii_list:        i-index of WWTP (only saved if cell is in model domain)
        jj_list:        j-index of WWTP (only saved if cell is in model domain)
        distance_list:  distance from grid cell to original WWTP location [m]  
    '''
    ii_list = []
    jj_list = []
    distance_list = []
    # if cell is in domain and it is a water cell (mask = 0),
    # then record information
    if cell_in_domain(I_ind,J_ind,len(X)-1,len(Y)-1):
        if mask_rho[I_ind,J_ind]:
            # record cell indices
            ii_list.append(I_ind)
            jj_list.append(J_ind)
            # record distance from center of cell to wwtp location
            xmeters, ymeters = zfun.ll2xy(X[I_ind], Y[J_ind], x, y)
            distance = np.sqrt(xmeters**2 + ymeters**2)
            distance_list.append(distance)
    return ii_list, jj_list, distance_list

def get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho):
    
    '''
    Function that check whether a grid cell is a coastal water cell.
    If it is coastal, returns updated row, col, idir, isign, and uv, 
    and the distance from the coastal cell to the river mouth. 

    River direction is decided based on the position of the closest 
    coastal land cell to the original river mouth lat/lon coordinates.

    Inputs:
        I_ind:      i-index of river mouth
        J_ind:      j-index of river mouth
        X:          domain lon coordinates
        Y:          domain lat coordinates
        x:          river mouth lon coordinates
        y:          river mouth lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        ii_list:        i-index of river (only saved if cell is a costal cell)
        jj_list:        j-index of river (only saved if cell is a costal cell)
        distance_list:  distance from grid cell to original river mouth location [m]  
        idir_list:      0 = E-W river, 1 = N-S river
        isign_list:     1 = flowing N/E, -1 = flowing S,W
        uv_list:        'u' = E-W river, 'v' = N-S river
    '''
    ii_list = []
    jj_list = []
    distance_list = []
    idir_list = []
    isign_list = []
    uv_list = []
    # if cell is in domain and it is a water cell (mask = 1),
    # then check if it is coastal
    if cell_in_domain(I_ind,J_ind,len(X)-1,len(Y)-1):
        if mask_rho[I_ind,J_ind]: # is water cell

            # Get coordinates of adjacent cells
            # north
            NI = I_ind
            NJ = J_ind + 1
            # south
            SI = I_ind
            SJ = J_ind - 1
            # west
            WI = I_ind - 1
            WJ = J_ind
            # east
            EI = I_ind + 1
            EJ = J_ind

            # tracker if there's a coastal cell or not
            is_coastal = False
            coastal_type = []

            # check if north cell is land
            if cell_in_domain(NI,NJ,len(X)-1,len(Y)-1):
                # check if north cell is land (meaning original cell is coastal)
                if not mask_rho[NI,NJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['N']
            # check if south cell is land
            if cell_in_domain(SI,SJ,len(X)-1,len(Y)-1):
                # check if south cell is land (meaning original cell is coastal)
                if not mask_rho[SI,SJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['S']
            # check if west cell is land
            if cell_in_domain(WI,WJ,len(X)-1,len(Y)-1):
                # check if west cell is land (meaning original cell is coastal)
                if not mask_rho[WI,WJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['W']
            # check if east cell is land
            if cell_in_domain(EI,EJ,len(X)-1,len(Y)-1):
                # check if east cell is land (meaning original cell is coastal)
                if not mask_rho[EI,EJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['E']

            # write data if the cell is coastal
            if is_coastal:
                # record distance from center of cell to source location
                xmeters, ymeters = zfun.ll2xy(X[I_ind], Y[J_ind], x, y)
                distance = np.sqrt(xmeters**2 + ymeters**2)

                # Figure out which adjacent land cell to use
                if len(coastal_type) == 1: # easy if there is only one land cell
                    source_side = coastal_type[0]
                # if there's more than one land cell, pick the one that is closest to the original source
                else:
                    dist2land = []
                    if 'N' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[NI], Y[NJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'S' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[SI], Y[SJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'W' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[WI], Y[WJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'E' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[EI], Y[EJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    # minimum distance
                    min_dist2land = np.min(dist2land)
                    # index of minimum distance
                    ind_dist2land = dist2land.index(min_dist2land)
                    # corresponding river direction
                    source_side = coastal_type[ind_dist2land]

                # get river direction values
                xoff = 0
                yoff = 0
                # offsets follow same convention as in pgrid code
                if source_side == 'N':
                        idir = 1
                        isign = -1
                        uv = 'v'
                elif source_side == 'S':
                        idir = 1
                        isign = 1
                        uv = 'v'
                        yoff = -1
                elif source_side == 'W':
                        idir = 0
                        isign = 1
                        uv = 'u'
                        xoff = -1
                elif source_side == 'E':
                        idir = 0
                        isign = -1
                        uv = 'u'

                # record distance
                distance_list.append(distance)
                # record cell indices
                ii_list.append(I_ind + xoff)
                jj_list.append(J_ind + yoff)
                # record cell flow directions
                idir_list.append(idir)
                isign_list.append(isign)
                uv_list.append(uv)

    return ii_list, jj_list, distance_list, idir_list, isign_list, uv_list

def get_nearest_coastal_cell_wwtp(sname,x,y,X,Y,mask_rho):
    """
    Algorithm that calculates the nearest coastal grid cell (row/col)
    to a point source,
    given its lat/lon coordinates, the lat/lon coordinates of the grid,
    and the mask of the rho-cells (to identify water vs land)

    Inputs:
        sname:      source name
        x:          WWTP lon coordinates
        y:          WWTP lat coordinates
        X:          domain lon coordinates
        Y:          domain lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        row:        row index of nearest coastal cell    
        col:        col index of nearest coastal cell
        ringnum:    ring that contains nearest coastal cell  
    """

    if in_domain(x, y, X, Y):

        # initialize a boolean to track whether a grid cell has been selected for point source
        ps_located = False

        # figure out in which grid cell the wwtp is located
        ix = zfun.find_nearest_ind(X,x)
        iy = zfun.find_nearest_ind(Y,y)
        wwtpII = np.array(ix)
        wwtpJJ = np.array(iy)

        # search for best grid cell to place wwtp-----------------------
    
        # First, check if the wwtp is located in the ocean
        if mask_rho[wwtpII,wwtpJJ]: # mask of 1 means water
            # point source located at same location as wwtp
            ringnum = 0
            psII = wwtpII
            psJJ = wwtpJJ
            ps_located = True

        # If on land, search in square rings around the source for the nearest coastal cell
        elif not mask_rho[wwtpII,wwtpJJ]:
            # initialize counters and arrays
            ringnum = 1
            II_list = []
            JJ_list = []
            dist_list = []
            min_dist = 1e9 # arbitratily large number so we can update with min value later

            while not ps_located:

                # record all water cells that are in top row of ring
                for iter in range(-ringnum,ringnum+1):
                    I_ind = wwtpII + iter
                    J_ind = wwtpJJ + ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in bottom row of ring
                for iter in range(-ringnum,ringnum+1):
                    I_ind = wwtpII + iter
                    J_ind = wwtpJJ - ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in left column of ring (exclude diagonals)
                for iter in range(-(ringnum-1),ringnum):
                    J_ind = wwtpJJ + iter
                    I_ind = wwtpII - ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in right column of ring (exclude diagonals)
                for iter in range(-(ringnum-1),ringnum):
                    J_ind = wwtpJJ + iter
                    I_ind = wwtpII + ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # calculate minimum distance coastal cell (if there were any coastal cells)
                if len(dist_list) > 0:
                    # minimum distance
                    min_dist_new = np.min(dist_list)
                    # index of minimum distance
                    min_dist_ind = dist_list.index(min_dist_new)
                    # corresponding i and j indices
                    min_II = II_list[min_dist_ind]
                    min_JJ = JJ_list[min_dist_ind]

                    # check if new min_dist is smaller than previous min_dist
                    if min_dist_new < min_dist:
                        # save new min dist and grid cell index
                        min_dist = min_dist_new
                        psII = min_II
                        psJJ = min_JJ
                        # check next ring if there is a closer coastal grid cell
                        ringnum +=1

                    # if new min_dist is larger than before, then we have found the nearest coastal grid cell
                    else:
                        ps_located = True

                # if all cells in ring are on land, then go into next ring
                else:
                    # iterate into next level ring
                    ringnum += 1

        row = psJJ
        col = psII
        return [row,col,ringnum]
    
    else:
        print(' >> excluding ' + sname)
        row = np.nan
        col = np.nan
        ringnum = np.nan
        return [row,col,ringnum]
        
def get_nearest_coastal_cell_riv(sname,x,y,X,Y,mask_rho):
    """
    Algorithm that calculates the nearest coastal grid cell (row/col)
    to a rivermouth,
    given its lat/lon coordinates, the lat/lon coordinates of the grid,
    and the mask of the rho-cells (to identify water vs land)
    
    Inputs:
        sname:      source name
        x:          river mouth lon coordinates
        y:          river mouth lat coordinates
        X:          domain lon coordinates
        Y:          domain lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        row:        row index of nearest coastal cell    
        col:        col index of nearest coastal cell
        idir:       0 = E-W river, 1 = N-S river
        isign:      1 = flowing N/E, -1 = flowing S,W
        uv:         'u' = E-W river, 'v' = N-S river
        ringnum:    ring that contains nearest coastal cell  
    """

    if in_domain(x, y, X, Y):

        # initialize a boolean to track whether a grid cell has been selected for point source
        rm_located = False

        # figure out in which grid cell the source is located
        ix = zfun.find_nearest_ind(X,x)
        iy = zfun.find_nearest_ind(Y,y)
        rivII = np.array(ix)
        rivJJ = np.array(iy)
    
        # search for best grid cell to place river mouth -----------------------
    
        # First, check if the river is already located in a coastal cell
        if mask_rho[rivII,rivJJ]: # mask of 1 means water
            ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(rivII,rivJJ,X,Y,x,y,mask_rho)
            # if list is not empty, then the cell is coastal
            if len(ii_list) > 0:
                # save the coastal information
                ringnum = 0
                rmII = ii_list[0]
                rmJJ = jj_list[0]
                idir = idir_list[0]
                isign = isign_list[0]
                uv = uv_list[0]
                rm_located = True

        # If original cell wasn't coastal,
        # Search in square rings around the source for the nearest coastal cell
        # initialize counters and arrays
        if not rm_located:
            ringnum = 1
            II_list = []
            JJ_list = []
            dist_list = []
            IDIR_list = []
            ISIGN_list = []
            UV_list = []
            min_dist = 1e9 # arbitratily large number so we can update with min value later

        while not rm_located:

            # record all water cells that are in top row of ring
            for iter in range(-ringnum,ringnum+1):
                I_ind = rivII + iter
                J_ind = rivJJ + ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in bottom row of ring
            for iter in range(-ringnum,ringnum+1):
                I_ind = rivII + iter
                J_ind = rivJJ - ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in left column of ring (exclude diagonals)
            for iter in range(-(ringnum-1),ringnum):
                J_ind = rivJJ + iter
                I_ind = rivII - ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in right column of ring (exclude diagonals)
            for iter in range(-(ringnum-1),ringnum):
                J_ind = rivJJ + iter
                I_ind = rivII + ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # calculate minimum distance coastal cell (if there were any cells in water)
            if len(dist_list) > 0:
                # minimum distance
                min_dist_new = np.min(dist_list)
                # index of minimum distance
                min_dist_ind = dist_list.index(min_dist_new)
                # corresponding i and j indices
                min_II = II_list[min_dist_ind]
                min_JJ = JJ_list[min_dist_ind]
                # corresponding river directions
                min_IDIR = IDIR_list[min_dist_ind]
                min_ISIGN = ISIGN_list[min_dist_ind]
                min_UV = UV_list[min_dist_ind]

                # check if new min_dist is smaller than previous min_dist
                if min_dist_new < min_dist:
                    # save new min dist and grid cell index
                    min_dist = min_dist_new
                    rmII = min_II
                    rmJJ = min_JJ
                    idir = min_IDIR
                    isign = min_ISIGN
                    uv = min_UV
                    # check next ring if there is a closer coastal grid cell
                    ringnum +=1

                # if new min_dist is larger than before, then we have found the nearest coastal grid cell
                else:
                    rm_located = True
                    # print('Originally on land')
                    # # print('New lat/lon:({},{})'.format(round(Y[psJJ],2),round(X[psII],2)))

            # if all cells in ring are on land, then go into next ring
            else:
                # iterate into next level ring
                ringnum += 1

        # save location of river mouth
        # print('New lat/lon:({},{})'.format(round(Y[rmJJ],2),round(X[rmII],2)))
        row = rmJJ
        col = rmII
        # print('Corresponding grid row/col:({},{})'.format(row,col))
        return [row,col,idir,isign,uv,ringnum]
    
    else:
        print(' >> excluding ' + sname)
        row = np.nan
        col = np.nan
        ringnum = np.nan
        idir = np.nan
        isign = np.nan
        uv = np.nan
        return [row,col,idir,isign,uv,ringnum]
