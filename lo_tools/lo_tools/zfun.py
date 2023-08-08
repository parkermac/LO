"""
This module contains utility functions for interpolation, filtering
and inspection.

Parker MacCready
"""

import numpy as np

def interp2(x, y, X, Y, U):
    """
    Interpolate field U(X,Y) to u(x,y).  All grids are required to be plaid and 2D
    """
    if is_plaid(x) and is_plaid(y) and is_plaid(X) and is_plaid(Y):
        # Much Faster than interp_scatttered_on_plaid()
        xi0, xi1, xf = get_interpolant(x[0,:], X[0,:])
        yi0, yi1, yf = get_interpolant(y[:,0], Y[:,0])
        NR, NC = x.shape
        XF, YF = np.meshgrid(xf, yf)
        # bi linear interpolation
        u00 = U[yi0,:][:,xi0]
        u10 = U[yi1,:][:,xi0]
        u01 = U[yi0,:][:,xi1]
        u11 = U[yi1,:][:,xi1]
        u = (1-YF)*((1-XF)*u00 + XF*u01) + YF*((1-XF)*u10 + XF*u11)
        return u
    else:
        print('grids not plaid')
        pass
        

def is_plaid(x):
    """
    Test if a numpy array is plaid.
    """
    if not isinstance(x, np.ndarray):
        return False
    elif not ((x[:,0]==x[:,1]).all() or (x[0,:]==x[1,:]).all()):
        return False
    else:
        return True

def interp_scattered_on_plaid(x, y, xvec, yvec, u):
    """
    Gets values of the field u at locations (x,y).

    All inputs and outputs are numpy arrays.

    Field u is defined on a plaid grid defined by vectors xvec and yvec.

    Locations (x,y) are defined by vectors whose elements specify
    individual points.  They are not a plaid grid, but can be
    scattered arbitrarily throughout the domain.

    Returns a vector ui the same length as x and y.

    """
    # get interpolants
    xi0, xi1, xf = get_interpolant(x,xvec)
    yi0, yi1, yf = get_interpolant(y,yvec)

    # bi linear interpolation
    u00 = u[yi0,xi0]
    u10 = u[yi1,xi0]
    u01 = u[yi0,xi1]
    u11 = u[yi1,xi1]
    ui = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)

    return ui

def get_interpolant(x, xvec, show_warnings=True):
    """
    Returns info to allow fast interpolation.

    Input:
    x = data position(s) [1-D numpy array]
    xvec = coordinate vector [1-D numpy array without nans]
        NOTE: xvec must be monotonically increasing
    
    Output: three 1-D numpy arrays of the same size as x
    i0 = index below [int]
    i1 = index above [int]
    fr = fraction [float]

    If the x is ON a point in xvec the default is to return
    the index of that point and the one above, with fr=0,
    unless it is the last point in which case it is the index
    of that point and the point below, with fr = 1.
    """
    
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    
    def itp_err(message='hi'):
        print('WARNING from get_interpolant(): ' + message)

    # input major error checking
    if isinstance(x, np.ndarray) and isinstance(xvec, np.ndarray):
        pass # input type is good
    else:
        itp_err('Inputs must be numpy arrays')
        return

    # some preconditioning of the input
    x = x.flatten()
    xvec = xvec.flatten()

    if show_warnings:
        # minor error checking
        if np.isnan(x).any():
            itp_err('nan found in x')
        if np.isnan(xvec).any():
            itp_err('nan found in xvec')
        if not np.all(np.diff(xvec) > 0):
            itp_err('xvec must be monotonic and increasing')

    nx = len(x)
    nxvec = len(xvec)

    X = x.reshape(nx, 1) # column vector
    xvec = xvec.reshape(1, nxvec)
    XVEC = xvec.repeat(nx, axis=0) # matrix

    # preallocate results arrays
    i0 = np.zeros(nx, dtype=int)
    i1 = np.zeros(nx, dtype=int)
    fr = np.zeros(nx, dtype=float)

    # calculate index columns
    mask = X >= XVEC
    # the above line broadcasts correctly even if nx = nxvec
    # because we forced X to be a column vector
    i0 = mask.sum(axis=1) - 1

    # these masks are used to handle values of x beyond the range of xvec
    lomask = i0 < 0
    himask = i0 > nxvec - 2
    i0[lomask] = 0
    i0[himask] = nxvec - 2
    i1 = i0 + 1

    # compute the fraction
    xvec0 = xvec[0,i0]
    xvec1 = xvec[0,i1]
    fr = (x - xvec0)/(xvec1 - xvec0)

    # fractions for out of range x
    fr[lomask] = np.nan
    fr[himask] = np.nan
    # override for the case where x = the last point of xvec
    fr[X[:,0]==XVEC[0,-1]] = 1.0

    return i0, i1, fr

def find_nearest(array, value):
    # gives the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return int(idx)

def filt_AB8d(data):
    """
    % 8/28/2013  Parker MacCready
    % ** use ONLY with hourly data! **
    %
    % This applies the Austin-Barth 8 day filter to a single vector.
    %
    % NOTE this may be different than their definition - it just returns a
    % weighted average of the data over the previous 8 days from time t,
    % with the weighting decaying from 1 to 1/e at t - 8 days.  There are 192
    % hours in 8 days.
    Input:
        assumed to be a 1D numpy array
    Output:
        a vector of the same size you started with,
        padded with NaN's at the ends.
    """
    fl = 8*24;
    filt = np.exp(np.linspace(-1,0,fl))
    filt = filt/filt.sum();
    smooth = np.nan*data
    for ii in range(fl+1, len(data)):
        smooth[ii] = (filt*data[ii-fl+1: ii + 1]).sum()
    return smooth
    
def lowpass(data, f='hanning', n=40, nanpad=True):
    """
    A replacement for almost all previous filter code.
    f = 'hanning' (default) or 'godin'
    
    Input: ND numpy array, any number of dimensions, with time on axis 0.
    
    Output: Array of the same size, filtered with Hanning window of length n,
        or the Godin filter (hourly data only) padded with nan's.
    """
    if n == 1:
        return data
    else:
        if f == 'hanning':
            filt = hanning_shape(n=n)
        elif f == 'godin':
            filt = godin_shape()
        else:
            print('ERROR in filt_general(): unsupported filter ' + f)
            filt = np.nan
        npad = np.floor(len(filt)/2).astype(int)
        sh = data.shape
        df = data.flatten('F')
        dfs = np.convolve(df, filt, mode = 'same')
        smooth = dfs.reshape(sh, order='F')
        # note that the indexing below defaults to being on axis 0,
        # and correctly broadcasts without having to mention the other axes
        if nanpad:
            smooth[:npad] = np.nan
            smooth[-npad:] = np.nan
        else:
            smooth[:npad] = data[:npad]
            smooth[-npad:] = data[-npad:]
        return smooth
    
def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def hanning_shape(n=40):
    """
    Returns a Hanning window of the specified length.
    """
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE

def ll2xy(lon, lat, lon0, lat0):
    """
    This converts lon, lat into meters relative to lon0, lat0.
    It should work for lon, lat scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y
    
def get_rc(NP):
    # figure out near-optimal numer of rows and columns for plotting
    NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
    NC = np.ceil(NP/NR).astype(int)
    return NR, NC
    
def get_irc(ii, NC):
    # get row and column of plot ii when there are NC columns
    ir = int(np.floor(ii/NC))
    ic = int(ii - NC*ir)
    return ir, ic
    
def fillit(a):
    # ensures a is an array with nan's for masked values
    # instead of a masked array
    if isinstance(a, np.ma.MaskedArray):
        a = a.filled(np.nan)
    return a
    
def boolean_string(s):
    # used by argparse
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

def dist(x, x1, y, y1):
    """
    Standard distance between two points, used by get_stairstep().
    """
    d = np.sqrt((x1-x)**2 + (y1-y)**2)
    return d
    
def dist_normal(x0, x1, y0, y1, xp, yp):
    """
    This finds the shortest distance from xp, yp to the segment
    defined by points (x0, y0) and (x1, y1).  Used by get_stairstep().
    """
    # make sure x is increasing (inside the function)
    if x1 < x0:
        (x0, x1) = (x1, x0)
        (y0, y1) = (y1, y0)
    dx = x1-x0
    dy = y1-y0
    # find xp2, yp2: the position of the point on the line
    # that is closest to xp, yp
    if dy == 0: # line is horizontal
        xp2 = xp
        yp2 = y0
    elif dx == 0: # line is vertical
        yp2 = yp
        xp2 = x0
    else: # The line is sloping: then we find the equation for
        # a line through xp, yp that is normal to the original line,
        # and then find the point xp2, yp2 where the two lines intersect.
        m = dy/dx
        b = y0 - m*x0
        bp = yp + (1/m)*xp
        xp2 = (bp - b)/(m + (1/m))
        yp2 = m*xp2 + b
    # and return the distance between the point xp, yp and the intersection point
    dxp = xp2 - xp
    dyp = yp2 - yp
    dist_n = np.sqrt(dxp**2 + dyp**2)
    return dist_n
    
def get_stairstep(x0, x1, y0, y1):
    """
    Brute force method of generating a "stairstep" path: always choosing
    the point with the shortest distance to the end will generate
    the straightest path.
    
    The inputs should be integers, representing indices into a grid.
    
    The outputs are arrays of indices that would be suitable for a river
    channel if we were working on the rho-grid on ROMS.
    """
    # check input
    if isinstance(x0,int) and isinstance(x1,int) and isinstance(y0,int) and isinstance(y1,int):
        pass
    else:
        print('Error from zfun.get_stairstep(): must pass ints.')
        return
    d = dist(x0, x1, y0, y1)
    xx = []
    yy = []
    xx.append(x0)
    yy.append(y0)
    x = x0
    y = y0
    while d > 0:
                
        # find distances of all 4 surrounding points that
        # are allowed choices to the end of line
        dn = dist(x, x1, y+1, y1)
        ds = dist(x, x1, y-1, y1)
        de = dist(x+1, x1, y, y1)
        dw = dist(x-1, x1, y, y1)
        # close_arr is an array that is positive for surrounding points
        # that are closer to the endpoint than the current point
        close_arr = np.array([d-dn, d-ds, d-de, d-dw])
        
        # find distances of all 4 surrounding points that
        # are allowed choices to closest point on line
        dn_n = dist_normal(x0, x1, y0, y1, x, y+1)
        ds_n = dist_normal(x0, x1, y0, y1, x, y-1)
        de_n = dist_normal(x0, x1, y0, y1, x+1, y)
        dw_n = dist_normal(x0, x1, y0, y1, x-1, y)
        # gather them in an array
        norm_arr = np.array([dn_n, ds_n, de_n, dw_n])
        
        # step_array is an array of the steps we take to make each
        # of the surrounding points
        step_arr = np.array([[0,1],[0,-1],[1,0],[-1,0]])
        
        # closer is a Boolean array that is True for all step choices
        # that got us closer to the target
        closer = close_arr >= 0
        
        # make shorter versions of step_arr and norm_arr that only include
        # choices that got us closer
        step_arr = step_arr[closer]
        norm_arr = norm_arr[closer]
        
        # then find the index of the remaining choice that was closest to the line
        imin = np.argmin(norm_arr)
        
        # take the step in the chosen direction
        step = step_arr[imin,:]
        x += step[0]
        y += step[1]
        xx.append(x)
        yy.append(y)
        
        # and update the distance so the loop knows when it is done
        d = dist(x, x1, y, y1)
        
    # pack results as arrays
    XX = np.array(xx)
    YY = np.array(yy)
    
    # check that result never steps more or less than one
    DD = np.abs(np.diff(XX)) + np.abs(np.diff(YY))
    if (DD==1).all():
        pass # the result passes this test
    else:
        print('Error in result stepsize')
    
    # check that endpoints are correct
    if (XX[0] != x0) or (XX[-1] != x1) or (YY[0] != y0) or (YY[-1] != y1):
        print('Error in result endpoints!')
        
    return XX, YY
    
def linefit(x,y):
    """
    Gives best fit line to vectors x and y, with statistics.
    """
    # do least-squares best fit to a line
    from scipy import stats
    BB = np.polyfit(x,y,1, full=True)[0]
    slope = BB[0] # slope
    y0 = BB[1] # y-intercept
    
    # fit line yfit(x)
    yfit = slope*x + y0
    
    # correlation coefficient
    r = np.corrcoef(x,y)[0,1]
    
    # Calculate the 95% confidence interval for the mean and trend
    # NOTE the call to the inverse cumulative Student's t-distribution:
    # stats.t.ppf(1-.025,100) which is identical to MATLAB tinv(1-.025,100).
    ci_pct = 95 # Confidence interval %
    ci_fac = (1 - ci_pct/100)/2 # = 0.025 for ci_pct = 95
    N = len(x)
    s_eps = np.sqrt(((y-yfit)**2).sum()/(N-2)) # standard error of the estimate
    s_x = x.std(ddof=1)
    s_y = y.std(ddof=1)
    # Condidence interval on mean: Emery and Thomson (3.8.6)
    ci_mean = s_y * stats.t.ppf(1-ci_fac,N-1) / np.sqrt(N)
    # Confidence interval on slope: Emery and Thomson (3.15.12b) w/typo corrected
    ci_trend = s_eps * stats.t.ppf(1-ci_fac,N-2) / np.sqrt((N-1)*s_x*s_x)
    
    return slope, y0, r, ci_mean, ci_trend
    







