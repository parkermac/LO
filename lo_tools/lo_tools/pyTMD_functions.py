"""
These functions are copied from the package https://github.com/tsutterley/pyTMD
which is mostly a direct python adaptation of the MATLAB tmd_toolbox.

The functions defined in this module are from a few different modules in pyTMD
(calc_astrol_longitudes.py, load_constituents.py, load_nodal_corrections.py, and
tidal_ellipse.py).  I gather them all in one module here just for convenience.

The use of these functions is to get the Greenwich phase for each tidal constituent,
relative to 1/1/1992, the reference time for the tpxo9 fields.  I also
use them to calculate the nodal corrections to phase and amplitude for a given
model day, and to calculate the ellipse parameters.

I have checked the results in a limited way:
* I checked to phase and amplitide numbers for M2 - they were identical.
* The calculated field patterns compared well with the cas6_v0/f2019.07.04/tide1 results.
The testing was done with LPM/tpxo/test1.py.
"""

import numpy as np

"""
load_constituent.py (04/2022)
Loads parameters for a given tidal constituent

CALLING SEQUENCE:
    amplitude,phase,omega,alpha,species = load_constituent(c)

INPUTS:
    c: tidal constituent ID

OUTPUT:
    amplitude: amplitude of equilibrium tide in m for tidal constituent
    phase: phase of tidal constituent
    omega: angular frequency of constituent in radians
    alpha: load love number of tidal constituent
    species: spherical harmonic dependence of quadrupole potential

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: add more constituents from OTPSnc and function docstrings
    Updated 09/2017: Rewritten in Python
"""

def load_constituent(c):
    """
    Loads parameters for a given tidal constituent

    Parameters
    ----------
    c: list
        tidal constituent ID

    Returns
    -------
    amplitude: float
        amplitude of equilibrium tide in m for tidal constituent
    phase: float
        phase of tidal constituent
    omega: float
        angular frequency of constituent in radians
    alpha: float
        load love number of tidal constituent
    species: float
        spherical harmonic dependence of quadrupole potential

    References
    ----------
    .. [1] Egbert and Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).
    """
    #-- constituents array that are included in tidal program
    cindex = ['m2','s2','k1','o1','n2','p1','k2','q1','2n2','mu2','nu2','l2',
        't2','j1','m1','oo1','rho1','mf','mm','ssa','m4','ms4','mn4','m6','m8',
        'mk3','s6','2sm2','2mk3']
    #-- species type (spherical harmonic dependence of quadrupole potential)
    species_all = np.array([2,2,1,1,2,1,2,1,2,2,2,2,2,1,1,1,1,0,0,0,0,0,0,
        0,0,0,0,0,0])
    #-- loading love number
    #-- alpha = correction factor for first order load tides
    alpha_all = np.array([0.693,0.693,0.736,0.695,0.693,0.706,0.693,0.695,0.693,
        0.693,0.693,0.693,0.693,0.695,0.695,0.695,0.695,0.693,0.693,0.693,
        0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693])
    #-- omega: angular frequency of constituent, in radians
    omega_all = np.array([1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,
        1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,1.352405e-04,
        1.355937e-04,1.382329e-04,1.431581e-04,1.452450e-04,7.556036e-05,
        7.028195e-05,7.824458e-05,6.531174e-05,0.053234e-04,0.026392e-04,
        0.003982e-04,2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,
        5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,2.081166e-04])
    #-- Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)
    #-- phases for each constituent are referred to the time when the phase of
    #-- the forcing for that constituent is zero on the Greenich meridian
    phase_all = np.array([1.731557546,0.000000000,0.173003674,1.558553872,
        6.050721243,6.110181633,3.487600001,5.877717569,4.086699633,
        3.463115091,5.427136701,0.553986502,0.052841931,2.137025284,
        2.436575100,1.929046130,5.254133027,1.756042456,1.964021610,
        3.487600001,3.463115091,1.731557546,1.499093481,5.194672637,
        6.926230184,1.904561220,0.000000000,4.551627762,3.809122439])
    #-- amplitudes of equilibrium tide in m
    # amplitude_all = np.array([0.242334,0.112743,0.141565,0.100661,0.046397,
    amplitude_all = np.array([0.2441,0.112743,0.141565,0.100661,0.046397,
        0.046848,0.030684,0.019273,0.006141,0.007408,0.008811,0.006931,0.006608,
        0.007915,0.007915,0.004338,0.003661,0.042041,0.022191,0.019567,0.,0.,0.,
        0.,0.,0.,0.,0.,0.])

    #-- map between input constituent and cindex
    j = [j for j,val in enumerate(cindex) if (val == c.lower())]
    #-- set the values for the constituent
    if j:
        amplitude, = amplitude_all[j]
        phase, = phase_all[j]
        omega, = omega_all[j]
        alpha, = alpha_all[j]
        species, = species_all[j]
    else:
        amplitude = 0.0; phase = 0.0; omega = 0.0; alpha = 0.0; species = 0
    #-- return the values for the constituent
    return (amplitude,phase,omega,alpha,species)
    
"""
calc_astrol_longitudes.py (04/2022)
Modification of ASTROL fortran subroutine by Richard Ray 03/1999

Computes the basic astronomical mean longitudes: s, h, p, N and PP
Note N is not N', i.e. N is decreasing with time.

Formulae for the period 1990--2010 were derived by David Cartwright
MEEUS and ASTRO5 formulae are from versions of Meeus's Astronomical Algorithms

CALLING SEQUENCE:
    s,h,p,N,PP = calc_astrol_longitudes(MJD, ASTRO5=True)

INPUTS:
    MJD: Modified Julian Day of input date

OUTPUTS:
    s: mean longitude of moon (degrees)
    h: mean longitude of sun (degrees)
    p: mean longitude of lunar perigee (degrees)
    N: mean longitude of ascending lunar node (degrees)
    PP: longitude of solar perigee (degrees)

OPTIONS:
    MEEUS: use additional coefficients from Meeus Astronomical Algorithms
    ASTRO5: use Meeus Astronomical coefficients as implemented in ASTRO5

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: change time variable names to not overwrite functions
    Updated 07/2020: added function docstrings
    Updated 07/2018: added option ASTRO5 to use coefficients from Richard Ray
        for use with the GSFC Global Ocean Tides (GOT) model
        added longitude of solar perigee (PP) as an additional output
    Updated 09/2017: added option MEEUS to use additional coefficients
        from Meeus Astronomical Algorithms to calculate mean longitudes
    Updated 09/2017: Rewritten in Python
    Rewritten in Matlab by Lana Erofeeva 2003
    Written by Richard Ray 12/1990
"""

#-- PURPOSE: calculate the sum of a polynomial function of time
def polynomial_sum(coefficients, t):
    """
    Calculates the sum of a polynomial function of time

    Parameters
    ----------
    coefficients: list
        leading coefficient of polynomials of increasing order
    t: float
        delta time in units for a given astronomical longitudes calculation
    """
    #-- convert time to array if importing a single value
    t = np.atleast_1d(t)
    return np.sum([c * (t ** i) for i,c in enumerate(coefficients)],axis=0)

#-- PURPOSE: compute the basic astronomical mean longitudes
def calc_astrol_longitudes(MJD, MEEUS=False, ASTRO5=False):
    """
    Computes the basic astronomical mean longitudes: s, h, p, N and PP

    Parameters
    ----------
    MJD: float
        Modified Julian Day (MJD) of input date
    MEEUS: bool, default False
        use additional coefficients from Meeus Astronomical Algorithms
    ASTRO5: bool, default False
        use Meeus Astronomical coefficients as implemented in ASTRO5

    Returns
    -------
    s: float
        mean longitude of moon (degrees)
    h: float
        mean longitude of sun (degrees)
    p: float
        mean longitude of lunar perigee (degrees)
    N: float
        mean longitude of ascending lunar node (degrees)
    PP: float
        longitude of solar perigee (degrees)
    """
    circle = 360.0
    if MEEUS:
        #-- convert from MJD to days relative to 2000-01-01T12:00:00
        T = MJD - 51544.5
        #-- mean longitude of moon
        lunar_longitude = np.array([218.3164591, 13.17639647754579,
            -9.9454632e-13, 3.8086292e-20, -8.6184958e-27])
        s = polynomial_sum(lunar_longitude,T)
        #-- mean longitude of sun
        solar_longitude = np.array([280.46645, 0.985647360164271,
            2.2727347e-13])
        h = polynomial_sum(solar_longitude,T)
        #-- mean longitude of lunar perigee
        lunar_perigee = np.array([83.3532430, 0.11140352391786447,
            -7.7385418e-12, -2.5636086e-19, 2.95738836e-26])
        p = polynomial_sum(lunar_perigee,T)
        #-- mean longitude of ascending lunar node
        lunar_node = np.array([125.0445550, -0.052953762762491446,
            1.55628359e-12, 4.390675353e-20, -9.26940435e-27])
        N = polynomial_sum(lunar_node,T)
        #-- mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T
    elif ASTRO5:
        #-- convert from MJD to centuries relative to 2000-01-01T12:00:00
        T = (MJD - 51544.5)/36525.0
        #-- mean longitude of moon (p. 338)
        lunar_longitude = np.array([218.3164477, 481267.88123421, -1.5786e-3,
             1.855835e-6, -1.53388e-8])
        s = polynomial_sum(lunar_longitude,T)
        #-- mean longitude of sun (p. 338)
        lunar_elongation = np.array([297.8501921, 445267.1114034, -1.8819e-3,
             1.83195e-6, -8.8445e-9])
        h = polynomial_sum(lunar_longitude-lunar_elongation,T)
        #-- mean longitude of lunar perigee (p. 343)
        lunar_perigee = np.array([83.3532465, 4069.0137287, -1.032e-2,
            -1.249172e-5])
        p = polynomial_sum(lunar_perigee,T)
        #-- mean longitude of ascending lunar node (p. 144)
        lunar_node = np.array([125.04452, -1934.136261, 2.0708e-3, 2.22222e-6])
        N = polynomial_sum(lunar_node,T)
        #-- mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T
    else:
        #-- convert from MJD to days relative to 2000-01-01T12:00:00
        #-- convert from Universal Time to Dynamic Time at 2000-01-01
        T = MJD - 51544.4993
        #-- mean longitude of moon
        s = 218.3164 + 13.17639648 * T
        #-- mean longitude of sun
        h = 280.4661 + 0.98564736 * T
        #-- mean longitude of lunar perigee
        p =  83.3535 + 0.11140353 * T
        #-- mean longitude of ascending lunar node
        N = 125.0445 - 0.05295377 * T
        #-- solar perigee at epoch 2000
        PP = 282.8

    #-- take the modulus of each
    s = np.mod(s, circle)
    h = np.mod(h, circle)
    p = np.mod(p, circle)
    N = np.mod(N, circle)

    #-- return as tuple
    return (s, h, p, N, PP)
    
"""
load_nodal_corrections.py (04/2022)
Calculates the nodal corrections for tidal constituents
Modification of ARGUMENTS fortran subroutine by Richard Ray 03/1999

CALLING SEQUENCE:
    pu,pf,G = load_nodal_corrections(MJD,constituents)

INPUTS:
    MJD: Modified Julian Day of input date
    constituents: tidal constituent IDs

OUTPUTS:
    pu,pf: nodal corrections for the constituents
    G: phase correction in degrees

OPTIONS:
    DELTAT: time correction for converting to Ephemeris Time (days)
    CORRECTIONS: use nodal corrections from OTIS/ATLAS or GOT models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes

REFERENCES:
    A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model
        time series", Advances in Water Resources, 12, (1989).
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2020: fix k1 for FES models
    Updated 08/2020: change time variable names to not overwrite functions
        update nodal corrections for FES models
    Updated 07/2020: added function docstrings.  add shallow water constituents
    Updated 09/2019: added netcdf option to CORRECTIONS option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
    Rewritten in Matlab by Lana Erofeeva 01/2003
    Written by Richard Ray 03/1999
"""

def load_nodal_corrections(MJD,constituents,DELTAT=0.0,CORRECTIONS='OTIS'):
    """
    Calculates the nodal corrections for tidal constituents

    Parameters
    ----------
    MJD: float
        modified julian day of input date
    constituents: list
        tidal constituent IDs
    DELTAT: float, default 0.0
        time correction for converting to Ephemeris Time (days)
    CORRECTIONS: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    PM EDIT: phase and amplitude were reversed in the original - a typo I assume
    pu: float
        nodal correction for the constituent phase [radians]
    pf: float
        nodal correction for the constituent amplitude [dimensionless]
    G: float
        phase correction in degrees
        PM NOTE: I don't know what this is for.  

    References
    ----------
    .. [1] Doodson and Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    .. [2] Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    .. [3] Foreman and Henry, "The harmonic analysis of tidal model time
        series", Advances in Water Resources, 12, (1989).
    .. [4] Egbert and Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).
    """

    #-- constituents array (not all are included in tidal program)
    cindex = ['sa','ssa','mm','msf','mf','mt','alpha1','2q1','sigma1','q1',
        'rho1','o1','tau1','m1','chi1','pi1','p1','s1','k1','psi1','phi1',
        'theta1','j1','oo1','2n2','mu2','n2','nu2','m2a','m2','m2b','lambda2',
        'l2','t2','s2','r2','k2','eta2','mns2','2sm2','m3','mk3','s3','mn4',
        'm4','ms4','mk4','s4','s5','m6','s6','s7','s8','m8','mks2','msqm','mtm',
        'n4','eps2','z0']

    #-- degrees to radians
    dtr = np.pi/180.0

    #-- set function for astronomical longitudes
    ASTRO5 = True if CORRECTIONS in ('GOT','FES') else False
    #-- convert from Modified Julian Dates into Ephemeris Time
    s,h,p,omega,pp = calc_astrol_longitudes(MJD+DELTAT, ASTRO5=ASTRO5)
    hour = (MJD % 1)*24.0
    t1 = 15.0*hour
    t2 = 30.0*hour
    nt = len(np.atleast_1d(MJD))

    #-- Determine equilibrium arguments
    arg = np.zeros((nt,60))
    arg[:,0] = h - pp #-- Sa
    arg[:,1] = 2.0*h #-- Ssa
    arg[:,2] = s - p #-- Mm
    arg[:,3] = 2.0*s - 2.0*h #-- MSf
    arg[:,4] = 2.0*s #-- Mf
    arg[:,5] = 3.0*s - p #-- Mt
    arg[:,6] = t1 - 5.0*s + 3.0*h + p - 90.0 #-- alpha1
    arg[:,7] = t1 - 4.0*s + h + 2.0*p - 90.0 #-- 2Q1
    arg[:,8] = t1 - 4.0*s + 3.0*h - 90.0 #-- sigma1
    arg[:,9] = t1 - 3.0*s + h + p - 90.0 #-- q1
    arg[:,10] = t1 - 3.0*s + 3.0*h - p - 90.0 #-- rho1
    arg[:,11] = t1 - 2.0*s + h - 90.0 #-- o1
    arg[:,12] = t1 - 2.0*s + 3.0*h + 90.0 #-- tau1
    arg[:,13] = t1 - s + h + 90.0 #-- M1
    arg[:,14] = t1 - s + 3.0*h - p + 90.0 #-- chi1
    arg[:,15] = t1 - 2.0*h + pp - 90.0 #-- pi1
    arg[:,16] = t1 - h - 90.0 #-- p1
    if CORRECTIONS in ('OTIS','ATLAS','netcdf'):
        arg[:,17] = t1 + 90.0 #-- s1
    elif CORRECTIONS in ('GOT','FES'):
        arg[:,17] = t1 + 180.0 #-- s1 (Doodson's phase)
    arg[:,18] = t1 + h + 90.0 #-- k1
    arg[:,19] = t1 + 2.0*h - pp + 90.0 #-- psi1
    arg[:,20] = t1 + 3.0*h + 90.0 #-- phi1
    arg[:,21] = t1 + s - h + p + 90.0 #-- theta1
    arg[:,22] = t1 + s + h - p + 90.0 #-- J1
    arg[:,23] = t1 + 2.0*s + h + 90.0 #-- OO1
    arg[:,24] = t2 - 4.0*s + 2.0*h + 2.0*p #-- 2N2
    arg[:,25] = t2 - 4.0*s + 4.0*h #-- mu2
    arg[:,26] = t2 - 3.0*s + 2.0*h + p #-- n2
    arg[:,27] = t2 - 3.0*s + 4.0*h - p #-- nu2
    arg[:,28] = t2 - 2.0*s + h + pp #-- M2a
    arg[:,29] = t2 - 2.0*s + 2.0*h #-- M2
    arg[:,30] = t2 - 2.0*s + 3.0*h - pp #-- M2b
    arg[:,31] = t2 - s + p + 180.0 #-- lambda2
    arg[:,32] = t2 - s + 2.0*h - p + 180.0 #-- L2
    arg[:,33] = t2 - h + pp #-- T2
    arg[:,34] = t2 #-- S2
    arg[:,35] = t2 + h - pp + 180.0 #-- R2
    arg[:,36] = t2 + 2.0*h #-- K2
    arg[:,37] = t2 + s + 2.0*h - pp #-- eta2
    arg[:,38] = t2 - 5.0*s + 4.0*h + p #-- MNS2
    arg[:,39] = t2 + 2.0*s - 2.0*h #-- 2SM2
    arg[:,40] = 1.5*arg[:,29] #-- M3
    arg[:,41] = arg[:,18] + arg[:,29] #-- MK3
    arg[:,42] = 3.0*t1 #-- S3
    arg[:,43] = arg[:,26] + arg[:,29] #-- MN4
    arg[:,44] = 2.0*arg[:,29] #-- M4
    arg[:,45] = arg[:,29] + arg[:,34] #-- MS4
    arg[:,46] = arg[:,29] + arg[:,36] #-- MK4
    arg[:,47] = 4.0*t1 #-- S4
    arg[:,48] = 5.0*t1 #-- S5
    arg[:,49] = 3.0*arg[:,29] #-- M6
    arg[:,50] = 3.0*t2 #-- S6
    arg[:,51] = 7.0*t1 #-- S7
    arg[:,52] = 4.0*t2 #-- S8
    #-- shallow water constituents
    arg[:,53] = 4.0*arg[:,29] #-- m8
    arg[:,54] = arg[:,29] + arg[:,36] - arg[:,34] #-- mks2
    arg[:,55] = 4.0*s - 2.0*h #-- msqm
    arg[:,56] = 3.0*s - p #-- mtm
    arg[:,57] = 2.0*arg[:,26] #-- n4
    arg[:,58] = t2 - 5.0*s + 4.0*h + p #-- eps2
    #-- mean sea level
    arg[:,59] = 0.0 #-- Z0

    #-- determine nodal corrections f and u
    sinn = np.sin(omega*dtr)
    cosn = np.cos(omega*dtr)
    sin2n = np.sin(2.0*omega*dtr)
    cos2n = np.cos(2.0*omega*dtr)
    sin3n = np.sin(3.0*omega*dtr)

    #-- set nodal corrections
    f = np.zeros((nt,60))
    u = np.zeros((nt,60))
    #-- determine nodal corrections f and u for each model type
    if CORRECTIONS in ('OTIS','ATLAS','netcdf'):
        f[:,0] = 1.0 #-- Sa
        f[:,1] = 1.0 #-- Ssa
        f[:,2] = 1.0 - 0.130*cosn #-- Mm
        f[:,3] = 1.0 #-- MSf
        f[:,4] = 1.043 + 0.414*cosn #-- Mf
        temp1 = (1.0 + 0.203*cosn + 0.040*cos2n)**2
        temp2 = (0.203*sinn + 0.040*sin2n)**2
        f[:,5] = np.sqrt(temp1 + temp2) #-- Mt
        f[:,6] = 1.0 #-- alpha1
        f[:,7] = np.sqrt((1.0 + 0.188*cosn)**2 + (0.188*sinn)**2) #-- 2Q1
        f[:,8] = f[:,7] #-- sigma1
        f[:,9] = f[:,7] #-- q1
        f[:,10] = f[:,7] #-- rho1
        temp1 = (1.0 + 0.189*cosn - 0.0058*cos2n)**2
        temp2 = (0.189*sinn - 0.0058*sin2n)**2
        f[:,11] = np.sqrt(temp1 + temp2) #-- O1
        f[:,12] = 1.0 #-- tau1
        #-- Doodson's
        # Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-omega)*dtr)
        # Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-omega)*dtr)
        #-- Ray's
        Mtmp1 = 1.36*np.cos(p*dtr) + 0.267*np.cos((p-omega)*dtr)
        Mtmp2 = 0.64*np.sin(p*dtr) + 0.135*np.sin((p-omega)*dtr)
        f[:,13] = np.sqrt(Mtmp1**2 + Mtmp2**2) #-- M1
        f[:,14] = np.sqrt((1.0+0.221*cosn)**2+(0.221*sinn)**2) #-- chi1
        f[:,15] = 1.0 #-- pi1
        f[:,16] = 1.0 #-- P1
        f[:,17] = 1.0 #-- S1
        temp1 = (1.0 + 0.1158*cosn - 0.0029*cos2n)**2
        temp2 = (0.1554*sinn - 0.0029*sin2n)**2
        f[:,18] = np.sqrt(temp1 + temp2) #-- K1
        f[:,19] = 1.0 #-- psi1
        f[:,20] = 1.0 #-- phi1
        f[:,21] = 1.0 #-- theta1
        f[:,22] = np.sqrt((1.0+0.169*cosn)**2 + (0.227*sinn)**2) #-- J1
        temp1 = (1.0 + 0.640*cosn + 0.134*cos2n)**2
        temp2 = (0.640*sinn + 0.134*sin2n)**2
        f[:,23] = np.sqrt(temp1 + temp2) #-- OO1
        temp1 = (1.0 - 0.03731*cosn + 0.00052*cos2n)**2
        temp2 = (0.03731*sinn - 0.00052*sin2n)**2
        f[:,24] = np.sqrt(temp1 + temp2) #-- 2N2
        f[:,25] = f[:,24] #-- mu2
        f[:,26] = f[:,24] #-- N2
        f[:,27] = f[:,24] #-- nu2
        f[:,28] = 1.0 #-- M2a
        f[:,29] = f[:,24] #-- M2
        f[:,30] = 1.0 #-- M2b
        f[:,31] = 1.0 #-- lambda2
        Ltmp1 = 1.0 - 0.25*np.cos(2*p*dtr) - 0.11*np.cos((2.0*p-omega)*dtr) - 0.04*cosn
        Ltmp2 = 0.25*np.sin(2*p*dtr) + 0.11*np.sin((2.0*p-omega)*dtr) + 0.04*sinn
        f[:,32] = np.sqrt(Ltmp1**2 + Ltmp2**2) #-- L2
        f[:,33] = 1.0 #-- T2
        f[:,34] = 1.0 #-- S2
        f[:,35] = 1.0 #-- R2
        temp1 = (1.0 + 0.2852*cosn + 0.0324*cos2n)**2
        temp2 = (0.3108*sinn + 0.0324*sin2n)**2
        f[:,36] = np.sqrt(temp1 + temp2) #-- K2
        f[:,37] = np.sqrt((1.0 + 0.436*cosn)**2 + (0.436*sinn)**2) #-- eta2
        f[:,38] = f[:,29]**2 #-- MNS2
        f[:,39] = f[:,29] #-- 2SM2
        f[:,40] = 1.0 #-- M3 (wrong)
        f[:,41] = f[:,18]*f[:,29] #-- MK3
        f[:,42] = 1.0 #-- S3
        f[:,43] = f[:,29]**2 #-- MN4
        f[:,44] = f[:,43] #-- M4
        f[:,45] = f[:,43] #-- MS4
        f[:,46] = f[:,29]*f[:,36] #-- MK4
        f[:,47] = 1.0 #-- S4
        f[:,48] = 1.0 #-- S5
        f[:,49] = f[:,29]**3 #-- M6
        f[:,50] = 1.0 #-- S6
        f[:,51] = 1.0 #-- S7
        f[:,52] = 1.0 #-- S8
        #-- shallow water constituents
        f[:,53] = f[:,29]**4 #-- m8
        f[:,54] = f[:,29]*f[:,36] #-- mks2
        f[:,55] = f[:,4] #-- msqm
        f[:,56] = f[:,4] #-- mtm
        f[:,57] = f[:,29]**2 #-- n4
        f[:,58] = f[:,29] #-- eps2
        #-- mean sea level
        f[:,59] = 1.0 #-- Z0

        u[:,0] = 0.0 #-- Sa
        u[:,1] = 0.0 #-- Ssa
        u[:,2] = 0.0 #-- Mm
        u[:,3] = 0.0 #-- MSf
        u[:,4] = -23.7*sinn + 2.7*sin2n - 0.4*sin3n #-- Mf
        temp1 = -(0.203*sinn + 0.040*sin2n)
        temp2 = (1.0 + 0.203*cosn + 0.040*cos2n)
        u[:,5] = np.arctan(temp1/temp2)/dtr #-- Mt
        u[:,6] = 0.0 #-- alpha1
        u[:,7] = np.arctan(0.189*sinn/(1.0 + 0.189*cosn))/dtr #-- 2Q1
        u[:,8] = u[:,7] #-- sigma1
        u[:,9] = u[:,7] #-- q1
        u[:,10] = u[:,7] #-- rho1
        u[:,11] = 10.8*sinn - 1.3*sin2n + 0.2*sin3n #-- O1
        u[:,12] = 0.0 #-- tau1
        u[:,13] = np.arctan2(Mtmp2,Mtmp1)/dtr #-- M1
        u[:,14] = np.arctan(-0.221*sinn/(1.0+0.221*cosn))/dtr #-- chi1
        u[:,15] = 0.0 #-- pi1
        u[:,16] = 0.0 #-- P1
        u[:,17] = 0.0 #-- S1
        temp1 = (-0.1554*sinn + 0.0029*sin2n)
        temp2 = (1.0 + 0.1158*cosn - 0.0029*cos2n)
        u[:,18] = np.arctan(temp1/temp2)/dtr #-- K1
        u[:,19] = 0.0 #-- psi1
        u[:,20] = 0.0 #-- phi1
        u[:,21] = 0.0 #-- theta1
        u[:,22] = np.arctan(-0.227*sinn/(1.0+0.169*cosn))/dtr #-- J1
        temp1 = -(0.640*sinn + 0.134*sin2n)
        temp2 = (1.0 + 0.640*cosn + 0.134*cos2n)
        u[:,23] = np.arctan(temp1/temp2)/dtr #-- OO1
        temp1 = (-0.03731*sinn + 0.00052*sin2n)
        temp2 = (1.0 - 0.03731*cosn + 0.00052*cos2n)
        u[:,24] = np.arctan(temp1/temp2)/dtr #-- 2N2
        u[:,25] = u[:,24] #-- mu2
        u[:,26] = u[:,24] #-- N2
        u[:,27] = u[:,24] #-- nu2
        u[:,28] = 0.0 #-- M2a
        u[:,29] = u[:,24] #-- M2
        u[:,30] = 0.0 #-- M2b
        u[:,31] = 0.0 #-- lambda2
        u[:,32] = np.arctan(-Ltmp2/Ltmp1)/dtr #-- L2
        u[:,33] = 0.0 #-- T2
        u[:,34] = 0.0 #-- S2
        u[:,35] = 0.0 #-- R2
        temp1 = -(0.3108*sinn+0.0324*sin2n)
        temp2 = (1.0 + 0.2852*cosn + 0.0324*cos2n)
        u[:,36] = np.arctan(temp1/temp2)/dtr #-- K2
        u[:,37] = np.arctan(-0.436*sinn/(1.0 + 0.436*cosn))/dtr #-- eta2
        u[:,38] = u[:,29]*2.0 #-- MNS2
        u[:,39] = u[:,29] #-- 2SM2
        u[:,40] = 1.50*u[:,29] #-- M3
        u[:,41] = u[:,29] + u[:,18] #-- MK3
        u[:,42] = 0.0 #-- S3
        u[:,43] = 2.0*u[:,29] #-- MN4
        u[:,44] = u[:,43] #-- M4
        u[:,45] = u[:,29] #-- MS4
        u[:,46] = u[:,29] + u[:,36] #-- MK4
        u[:,47] = 0.0 #-- S4
        u[:,48] = 0.0 #-- S5
        u[:,49] = 3.0*u[:,29] #-- M6
        u[:,50] = 0.0 #-- S6
        u[:,51] = 0.0 #-- S7
        u[:,52] = 0.0 #-- S8
        #-- mean sea level
        u[:,59] = 0.0 #-- Z0

    elif CORRECTIONS in ('FES',):
        #-- additional astronomical terms for FES models
        II = np.arccos(0.913694997 - 0.035692561*np.cos(omega*dtr))
        at1 = np.arctan(1.01883*np.tan(omega*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(omega*dtr/2.0))
        xi = -at1 - at2 + omega*dtr
        xi[xi > np.pi] -= 2.0*np.pi
        nu = at1 - at2
        I2 = np.tan(II/2.0)
        Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(p - xi)) + 36.0*(I2**4))
        P2 = np.sin(2.0*(p - xi))
        Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(p - xi))
        R = np.arctan(P2/Q2)
        P_prime = np.sin(2.0*II)*np.sin(nu)
        Q_prime = np.sin(2.0*II)*np.cos(nu) + 0.3347
        nu_prime = np.arctan(P_prime/Q_prime)
        P_sec = (np.sin(II)**2)*np.sin(2.0*nu)
        Q_sec = (np.sin(II)**2)*np.cos(2.0*nu) + 0.0727
        nu_sec = 0.5*np.arctan(P_sec/Q_sec)

        f[:,0] = 1.0 #-- Sa
        f[:,1] = 1.0 #-- Ssa
        f[:,2] = (2.0/3.0 - np.power(np.sin(II),2.0))/0.5021 #-- Mm
        f[:,3] = 1.0 #-- MSf
        f[:,4] = np.power(np.sin(II),2.0)/0.1578  #-- Mf
        f[:,7] = np.sin(II)*(np.cos(II/2.0)**2)/0.38 #-- 2Q1
        f[:,8] = f[:,7] #-- sigma1
        f[:,9] = f[:,7] #-- q1
        f[:,10] = f[:,7] #-- rho1
        f[:,11] = f[:,7] #-- O1
        #-- Ray's
        Mtmp1 = 1.36*np.cos(p*dtr) + 0.267*np.cos((p-omega)*dtr)
        Mtmp2 = 0.64*np.sin(p*dtr) + 0.135*np.sin((p-omega)*dtr)
        f[:,13] = np.sqrt(Mtmp1**2 + Mtmp2**2) #-- M1
        f[:,14] = np.sin(2.0*II) / 0.7214 #-- chi1
        f[:,15] = 1.0 #-- pi1
        f[:,16] = 1.0 #-- P1
        f[:,17] = 1.0 #-- S1
        temp1 = 0.8965*np.power(np.sin(2.0*II),2.0)
        temp2 = 0.6001*np.sin(2.0*II)*np.cos(nu)
        f[:,18] = np.sqrt(temp1 + temp2 + 0.1006) #-- K1
        f[:,19] = 1.0 #-- psi1
        f[:,20] = 1.0 #-- phi1
        f[:,21] = f[:,14] #-- theta1
        f[:,22] = f[:,14] #-- J1
        f[:,23] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 #-- OO1
        f[:,24] = np.power(np.cos(II/2.0),4.0)/0.9154 #-- 2N2
        f[:,25] = f[:,24] #-- mu2
        f[:,26] = f[:,24] #-- N2
        f[:,27] = f[:,24] #-- nu2
        f[:,28] = 1.0 #-- M2a
        f[:,29] = f[:,24] #-- M2
        f[:,30] = 1.0 #-- M2b
        f[:,31] = f[:,29] #-- lambda2
        f[:,32] = f[:,29]*Ra1 #-- L2
        f[:,33] = 1.0 #-- T2
        f[:,34] = 1.0 #-- S2
        f[:,35] = 1.0 #-- R2
        temp1 = 19.0444 * np.power(np.sin(II),4.0)
        temp2 = 2.7702 * np.power(np.sin(II),2.0) * np.cos(2.0*nu)
        f[:,36] = np.sqrt(temp1 + temp2 + 0.0981) #-- K2
        f[:,37] = np.power(np.sin(II),2.0)/0.1565 #-- eta2
        f[:,38] = f[:,29]**2 #-- MNS2
        f[:,39] = f[:,29] #-- 2SM2
        f[:,40] = np.power(np.cos(II/2.0), 6.0) / 0.8758 #-- M3
        f[:,41] = f[:,18]*f[:,29] #-- MK3
        f[:,42] = 1.0 #-- S3
        f[:,43] = f[:,29]**2 #-- MN4
        f[:,44] = f[:,43] #-- M4
        f[:,45] = f[:,29] #-- MS4
        f[:,46] = f[:,29]*f[:,36] #-- MK4
        f[:,47] = 1.0 #-- S4
        f[:,48] = 1.0 #-- S5
        f[:,49] = f[:,29]**3 #-- M6
        f[:,50] = 1.0 #-- S6
        f[:,51] = 1.0 #-- S7
        f[:,52] = 1.0 #-- S8
        #-- shallow water constituents
        f[:,53] = f[:,29]**4 #-- m8
        f[:,54] = f[:,29]*f[:,36] #-- mks2
        f[:,55] = f[:,4] #-- msqm
        f[:,56] = f[:,4] #-- mtm
        f[:,57] = f[:,29]**2 #-- n4
        f[:,58] = f[:,29] #-- eps2
        #-- mean sea level
        f[:,59] = 1.0 #-- Z0

        u[:,0] = 0.0 #-- Sa
        u[:,1] = 0.0 #-- Ssa
        u[:,2] = 0.0 #-- Mm
        u[:,3] = (2.0*xi - 2.0*nu)/dtr #-- MSf
        u[:,4] = -2.0*xi/dtr #-- Mf
        u[:,7] = (2.0*xi - nu)/dtr #-- 2Q1
        u[:,8] = u[:,7] #-- sigma1
        u[:,9] = u[:,7] #-- q1
        u[:,10] = u[:,7] #-- rho1
        u[:,11] = u[:,7] #-- O1
        u[:,13] = np.arctan2(Mtmp2,Mtmp1)/dtr #-- M1
        u[:,14] = -nu/dtr #-- chi1
        u[:,15] = 0.0 #-- pi1
        u[:,16] = 0.0 #-- P1
        u[:,17] = 0.0 #-- S1
        u[:,18] = -nu_prime/dtr #-- K1
        u[:,19] = 0.0 #-- psi1
        u[:,20] = 0.0 #-- phi1
        u[:,21] = -nu/dtr #-- theta1
        u[:,22] = u[:,21] #-- J1
        u[:,23] = (-2.0*xi - nu)/dtr #-- OO1
        u[:,24] = (2.0*xi - 2.0*nu)/dtr #-- 2N2
        u[:,25] = u[:,24] #-- mu2
        u[:,26] = u[:,24] #-- N2
        u[:,27] = u[:,24] #-- nu2
        u[:,29] = u[:,24] #-- M2
        u[:,31] = (2.0*xi - 2.0*nu)/dtr #-- lambda2
        u[:,32] = (2.0*xi - 2.0*nu - R)/dtr #-- L2
        u[:,33] = 0.0 #-- T2
        u[:,34] = 0.0 #-- S2
        u[:,35] = 0.0 #-- R2
        u[:,36] = -2.0*nu_sec/dtr #-- K2
        u[:,37] = -2.0*nu/dtr #-- eta2
        u[:,38] = (4.0*xi - 4.0*nu)/dtr #-- mns2
        u[:,39] = (2.0*xi - 2.0*nu)/dtr #-- 2SM2
        u[:,40] = (3.0*xi - 3.0*nu)/dtr #-- M3
        u[:,41] = (2.0*xi - 2.0*nu - 2.0*nu_prime)/dtr #-- MK3
        u[:,42] = 0.0 #-- S3
        u[:,43] = (4.0*xi - 4.0*nu)/dtr #-- MN4
        u[:,44] = (4.0*xi - 4.0*nu)/dtr #-- M4
        u[:,45] = (2.0*xi - 2.0*nu)/dtr  #-- MS4
        u[:,46] = (2.0*xi - 2.0*nu - 2.0*nu_sec)/dtr #-- MK4
        u[:,47] = 0.0 #-- S4
        u[:,48] = 0.0 #-- S5
        u[:,49] = (6.0*xi - 6.0*nu)/dtr #-- M6
        u[:,50] = 0.0 #-- S6
        u[:,51] = 0.0 #-- S7
        u[:,52] = 0.0 #-- S8
        #-- shallow water constituents
        u[:,53] = (8.0*xi - 8.0*nu)/dtr #-- m8
        u[:,54] = (2.0*xi - 2.0*nu - 2.0*nu_sec)/dtr #-- mks2
        u[:,55] = u[:,4] #-- msqm
        u[:,56] = u[:,4] #-- mtm
        u[:,57] = (4.0*xi - 4.0*nu)/dtr #-- n4
        u[:,58] = u[:,29] #-- eps2
        #-- mean sea level
        u[:,59] = 0.0 #-- Z0

    elif CORRECTIONS in ('GOT',):
        f[:,9] = 1.009 + 0.187*cosn - 0.015*cos2n#-- Q1
        f[:,11] = f[:,9]#-- O1
        f[:,16] = 1.0 #-- P1
        f[:,17] = 1.0 #-- S1
        f[:,18] = 1.006 + 0.115*cosn - 0.009*cos2n#-- K1
        f[:,26] = 1.000 - 0.037*cosn#-- N2
        f[:,29] = f[:,26]#-- M2
        f[:,34] = 1.0 #-- S2
        f[:,36] = 1.024 + 0.286*cosn + 0.008*cos2n#-- K2
        f[:,44] = f[:,29]**2#-- M4

        u[:,9] = 10.8*sinn - 1.3*sin2n#-- Q1
        u[:,11] = u[:,9]#-- O1
        u[:,16] = 0.0 #-- P1
        u[:,17] = 0.0 #-- S1
        u[:,18] = -8.9*sinn + 0.7*sin2n#-- K1
        u[:,26] = -2.1*sinn#-- N2
        u[:,29] = u[:,26]#-- M2
        u[:,34] = 0.0 #-- S2
        u[:,36] = -17.7*sinn + 0.7*sin2n#-- K2
        u[:,44] = -4.2*sinn#-- M4

    #-- take pu,pf,G for the set of given constituents
    nc = len(constituents)
    pu = np.zeros((nt,nc))
    pf = np.zeros((nt,nc))
    G = np.zeros((nt,nc))
    for i,c in enumerate(constituents):
        #-- map between given constituents and supported in tidal program
        j, = [j for j,val in enumerate(cindex) if (val == c)]
        pu[:,i] = u[:,j]*dtr
        pf[:,i] = f[:,j]
        G[:,i] = arg[:,j]

    #-- return values as tuple
    return (pu,pf,G)
    
"""
tidal_ellipse.py (04/2022)
Expresses the amplitudes and phases for the u and v components in terms of
    four ellipse parameters using Foreman's formula

CALLING SEQUENCE:
    umajor,uminor,uincl,uphase = tidal_ellipse(u,v)

INPUTS:
    u: zonal current (EW)
    v: meridional current (NS)

OUTPUTS:
    umajor: amplitude of the semimajor semi-axis
    uminor: amplitude of the semiminor semi-axis
    uincl: angle of inclination of the northern semimajor semi-axis
    uphase: phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

REFERENCE:
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model time
        series", Advances in Water Resources, 12(3), 109-120, (1989).
        https://doi.org/10.1016/0309-1708(89)90017-1

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2020
"""

def tidal_ellipse(u,v):
    """
    Expresses the amplitudes and phases for the u and v components in terms of
    four ellipse parameters using Foreman's formula

    Parameters
    ----------
    u: float
        zonal current (EW)
    v: float
        meridional current (NS)

    Returns
    -------
    umajor: float
        amplitude of the semimajor semi-axis
    uminor: float
        amplitude of the semiminor semi-axis
    uincl: float
        angle of inclination of the northern semimajor semi-axis
    uphase: float
        phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

    References
    ----------
    .. [1] Foreman and Henry, "The harmonic analysis of tidal model time
        series", Advances in Water Resources, 12(3), 109-120, (1989).
        https://doi.org/10.1016/0309-1708(89)90017-1

    """
    #-- change to polar coordinates
    t1p = u.real - v.imag
    t2p = v.real + u.imag
    t1m = u.real + v.imag
    t2m = v.real - u.imag
    #-- ap, am: amplitudes of positively and negatively rotating vectors
    ap = np.sqrt(t1p**2 + t2p**2)/2.0
    am = np.sqrt(t1m**2 + t2m**2)/2.0
    #-- ep, em: phases of positively and negatively rotating vectors
    ep = 180.0*np.arctan2(t2p, t1p)/np.pi
    ep[ep < 0.0] += 360.0
    em = 180.0*np.arctan2(t2m, t1m)/np.pi
    em[em < 0.0] += 360.0
    #-- determine the amplitudes of the semimajor and semiminor axes
    #-- using Foreman's formula
    umajor = (ap + am)
    uminor = (ap - am)
    #-- determine the inclination and phase using Foreman's formula
    uincl = 0.5 * (em + ep)
    uincl[uincl > 180.0] -= 180.0
    uphase = -0.5*(ep - em)
    uphase[uphase < 0.0] += 360.0
    uphase[uphase >= 360.0] -= 360.0
    #-- return values
    return (umajor, uminor, uincl, uphase)