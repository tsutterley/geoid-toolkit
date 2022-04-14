#!/usr/bin/env python
u"""
real_potential.py
Written by Tyler Sutterley (04/2022)
Calculates the real potential at a given latitude and height using
    coefficients from a gravity model

CALLING SEQUENCE:
    W, dW_dr, dW_dtheta = real_potential(lat, lon, h, clm, slm, lmax, R, GM)

INPUT:
    latitude: latitude in degrees
    longitude: longitude in degrees
    height: height above reference ellipsoid in meters
    clm: cosine spherical harmonics for a gravity model
    slm: sin spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree
    R: average radius used in gravity model
    GM: geocentric gravitational constant used in gravity model

OPTIONS:
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

OUTPUT:
    W: real potential at height h
    dW_dr: derivative of real potential with respect to radius
    dW_dtheta: derivative of real potential with respect to theta

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

REFERENCE:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4
    Barthelmes, "Definition of Functionals of the Geopotential and Their
        Calculation from Spherical Harmonic Models", STR09/02 (2009)
        http://icgem.gfz-potsdam.de/str-0902-revised.pdf
    Moazezi and Zomorrodian, "GGMCalc a software for calculation of the geoid
        undulation and the height anomaly using the iteration method, and
        classical gravity anomaly", Earth Science Informatics (2012)
        https://doi.org/10.1007/s12145-012-0102-2
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        https://doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
        changed dtypes to long double for high degree and order models
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
from geoid_toolkit.gauss_weights import gauss_weights

def real_potential(lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=0):
    """
    Calculates the real potential at a given latitude, longitude and
    height using coefficients from a gravity model following
    [Barthelmes2013]_ and [Moazezi2012]_

    Parameters
    ----------
    lat: float
        latitude in degrees
    lon: float
        longitude in degrees
    h: float
        ellipsoidal height in meters
    refell: str
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model
    clm: float
        cosine spherical harmonics for a gravity model
    slm: float
        sine spherical harmonics for a gravity model
    lmax: int
        maximum spherical harmonic degree
    R: float
        average radius used in gravity model
    GM: float
        geocentric gravitational constant used in gravity model
    GAUSS: float, default 0
        Gaussian Smoothing Radius in km

    Returns
    -------
    W: float
        real potential at height h
    dW_dr: float
        derivative of real potential with respect to radius
    dW_dtheta: float
        derivative of real potential with respect to theta

    References
    ----------
    .. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the
        Geopotential and Their Calculation from Spherical Harmonic Models",
        *GeoForschungsZentrum Scientific Technical Report*, STR09/02, (2013).
        `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_
    .. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz,
        *Physical Geodesy*, 2nd Edition, 403 pp., (2006).
        `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
    .. [Moazezi2012] S. Moazezi and H. Zomorrodian,
        "GGMCalc a software for calculation of the geoid undulation and the height
        anomaly using the iteration method, and classical gravity anomaly",
        *Earth Science Informatics*, 5, 123--136, (2012).
        `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_
    """

    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = np.longdouble(ellip['a'])
    ecc1 = np.longdouble(ellip['ecc1'])

    #-- convert from geodetic latitude to geocentric latitude
    latitude_geodetic_rad = (np.pi*lat/180.0).astype(np.longdouble)
    longitude_rad = (np.pi*lon/180.0).astype(np.longdouble)
    #-- prime vertical radius of curvature
    N = a/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(longitude_rad)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    #-- height of the observation point above the ellipsoid
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    latitude_geocentric_rad = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))
    #-- number of observations
    nlat = len(lat)
    #-- sin and cos of latitude
    t = np.sin(latitude_geocentric_rad)
    u = np.cos(latitude_geocentric_rad)
    #-- radius ratio
    q = (R / rr).astype(np.longdouble)

    #-- smooth the global gravity field with a Gaussian function
    if (GAUSS != 0):
        wt = 2.0*np.pi*gauss_weights(GAUSS,lmax)
        for l in range(0,lmax+1):
            clm[l,:] = clm[l,:]*wt[l]
            slm[l,:] = slm[l,:]*wt[l]

    #-- calculate clenshaw summations
    s_m_c = np.zeros((nlat,lmax*2+2),dtype=np.longdouble)
    ds_m_dr_c = np.zeros((nlat,lmax*2+2),dtype=np.longdouble)
    for m in range(lmax, -1, -1):
        s_m_c[:,2*m:2*m+2] = clenshaw_s_m(t,q,m,clm,slm,lmax)
        ds_m_dr_c[:,2*m:2*m+2] = clenshaw_ds_m_dr(t,q,m,clm,slm,lmax)

    #-- calculate cos phi
    cos_phi_2 = 2.0*np.cos(longitude_rad)
    #-- matrix of cos/sin m*phi (longitude_rad) summation
    cos_m_phi = np.zeros((nlat,lmax+2),dtype=np.longdouble)
    sin_m_phi = np.zeros((nlat,lmax+2),dtype=np.longdouble)
    #-- initialize matrix with values at lmax+1 and lmax
    cos_m_phi[:,lmax+1] = np.cos(np.longdouble(lmax + 1)*longitude_rad)
    sin_m_phi[:,lmax+1] = np.sin(np.longdouble(lmax + 1)*longitude_rad)
    cos_m_phi[:,lmax] = np.cos(np.longdouble(lmax)*longitude_rad)
    sin_m_phi[:,lmax] = np.sin(np.longdouble(lmax)*longitude_rad)
    #-- calculate summation
    s_m = s_m_c[:,2*lmax]*cos_m_phi[:,lmax] + s_m_c[:,2*lmax+1]*sin_m_phi[:,lmax]
    ds_m_dr = ds_m_dr_c[:,2*lmax]*cos_m_phi[:,lmax] + \
        ds_m_dr_c[:,2*lmax+1]*sin_m_phi[:,lmax]

    #-- iterate to calculate complete summation
    for m in range(lmax-1, 0, -1):
        cos_m_phi[:,m] = cos_phi_2*cos_m_phi[:,m+1] - cos_m_phi[:,m+2]
        sin_m_phi[:,m] = cos_phi_2*sin_m_phi[:,m+1] - sin_m_phi[:,m+2]
        a_m = np.sqrt((2.0*np.longdouble(m)+3.0)/(2.0*np.longdouble(m)+2.0))
        s_m = a_m*u*q*s_m + s_m_c[:,2*m]*cos_m_phi[:,m] + \
            s_m_c[:,2*m+1]*sin_m_phi[:,m]
        ds_m_dr = a_m*u*q*ds_m_dr + ds_m_dr_c[:,2*m]*cos_m_phi[:,m] + \
            ds_m_dr_c[:,2*m+1]*sin_m_phi[:,m]

    #-- add the final terms
    s_m = np.sqrt(3.0)*u*q*s_m + s_m_c[:,0]
    ds_m_dr = np.sqrt(3.0)*u*q*ds_m_dr + ds_m_dr_c[:,0]
    #-- compute the real potential and derivatives
    W = (GM/rr) * s_m
    dW_dr = (GM / (rr**2.0)) * ds_m_dr
    #-- return the real potential and derivatives
    return (W,dW_dr)

#-- PURPOSE: compute Clenshaw summation of the fully normalized associated
#-- Legendre's function for constant order m
def clenshaw_s_m(t, q, m, clm1, slm1, lmax):
    #-- allocate for output matrix
    N = len(t)
    s_m = np.zeros((N,2),dtype=np.longdouble)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.longdouble)
    slm = scalef*slm1.astype(np.longdouble)
    #-- convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if (m == lmax):
        s_m[:,0] = np.copy(clm[lmax,lmax])
        s_m[:,1] = np.copy(slm[lmax,lmax])
    elif (m == (lmax-1)):
        a_lm = t*q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        s_m[:,0] = a_lm*clm[lmax,lmax-1] + clm[lmax-1,lmax-1]
        s_m[:,1] = a_lm*slm[lmax,lmax-1] + slm[lmax-1,lmax-1]
    elif ((m <= (lmax-2)) and (m >= 1)):
        s_mm_c_pre_2 = np.copy(clm[lmax,m])
        s_mm_s_pre_2 = np.copy(slm[lmax,m])
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t*q
        s_mm_c_pre_1 = a_lm*s_mm_c_pre_2 + clm[lmax-1,m]
        s_mm_s_pre_1 = a_lm*s_mm_s_pre_2 + slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*t*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))*q**2
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,m]
            s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + slm[l,m]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_s_pre_2 = np.copy(s_mm_s_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
            s_mm_s_pre_1 = np.copy(s_mm_s)
        s_m[:,0] = np.copy(s_mm_c)
        s_m[:,1] = np.copy(s_mm_s)
    elif (m == 0):
        s_mm_c_pre_2 = np.copy(clm[lmax,0])
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))*t*q
        s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1)*(ll+1)))*t*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+1.)*(ll+1.))/((ll+2)*(ll+2)*(2.*ll+1.)))*q**2
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,0]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
        s_m[:,0] = np.copy(s_mm_c)
    #-- return s_m rescaled with scalef
    return s_m/scalef

#-- PURPOSE: compute Clenshaw summation of derivative with respect to latitude
#-- of the fully normalized associated Legendre's function for constant order m
def clenshaw_ds_m(t, u, q, m, clm1, slm1, lmax):
    #-- allocate for output matrix
    N = len(t)
    s_m = np.zeros((N,2),dtype=np.longdouble)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.longdouble)
    slm = scalef*slm1.astype(np.longdouble)
    #-- convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if (m == lmax):
        s_m[:,0] = mm*t*u*clm[lmax,lmax]
        s_m[:,1] = mm*t*u*slm[lmax,lmax]
    elif (m == (lmax-1)):
        a_lm = q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        s_dot_mm_c = a_lm * clm[lmax,lmax-1]
        s_dot_mm_s = a_lm * slm[lmax,lmax-1]
        s_mm_c = a_lm*t*clm[lmax,lmax-1] + clm[lmax-1,lmax-1]
        s_mm_s = a_lm*t*slm[lmax,lmax-1] + slm[lmax-1,lmax-1]
        s_m[:,0] = mm * t * u * s_mm_c - u * s_dot_mm_c
        s_m[:,1] = mm * t * u * s_mm_c - u * s_dot_mm_c
    elif ((m <= (lmax-2)) and (m >= 1)):
        s_mm_c_pre_2 = np.copy(clm[lmax,m])
        s_mm_s_pre_2 = np.copy(slm[lmax,m])
        s_dot_mm_c_pre_2 = 0.0
        s_dot_mm_s_pre_2 = 0.0
        a_lm = q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        s_dot_mm_c_pre_1 = a_lm * s_mm_c_pre_2
        s_dot_mm_s_pre_1 = a_lm * s_mm_s_pre_2
        s_mm_c_pre_1 = a_lm*t*s_mm_c_pre_2 + clm[lmax-1,m]
        s_mm_s_pre_1 = a_lm*t*s_mm_s_pre_2 + slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))*q**2
            s_dot_mm_c = a_lm * (s_dot_mm_c_pre_1 * t + s_mm_c_pre_1) - b_lm * s_dot_mm_c_pre_2
            s_dot_mm_s = a_lm * (s_dot_mm_s_pre_1 * t + s_mm_s_pre_1) - b_lm * s_dot_mm_s_pre_2
            s_mm_c = a_lm * t * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,m]
            s_mm_s = a_lm * t * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + slm[l,m]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_s_pre_2 = np.copy(s_mm_s_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
            s_mm_s_pre_1 = np.copy(s_mm_s)
            s_dot_mm_c_pre_2 = np.copy(s_dot_mm_c_pre_1)
            s_dot_mm_s_pre_2 = np.copy(s_dot_mm_s_pre_1)
            s_dot_mm_c_pre_1 = np.copy(s_dot_mm_c)
            s_dot_mm_s_pre_1 = np.copy(s_dot_mm_s)
        s_m[:,0] = mm * t * u * s_mm_c - u * s_dot_mm_c
        s_m[:,1] = mm * t * u * s_mm_s - u * s_dot_mm_s
    elif (m == 0):
        s_mm_c_pre_2 = np.copy(clm[lmax,0])
        a_lm = q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))
        s_dot_mm_c_pre_1 = a_lm * s_mm_c_pre_2
        s_mm_c_pre_1 = a_lm * t * s_mm_c_pre_2 + clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1)*(ll+1)))*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+1.)*(ll+1.))/((ll+2)*(ll+2)*(2.*ll+1.)))*q**2
            s_dot_mm_c = a_lm * (s_dot_mm_c_pre_1 * t + s_mm_c_pre_1) - b_lm * s_dot_mm_c_pre_2
            s_mm_c = a_lm * t * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,0]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
            s_dot_mm_c_pre_2 = np.copy(s_dot_mm_c_pre_1)
            s_dot_mm_c_pre_1 = np.copy(s_dot_mm_c)
        s_m[:,0] = - u * s_dot_mm_c
    #-- return s_m rescaled with scalef
    return s_m/scalef

#-- PURPOSE: compute Clenshaw summation of derivative with respect to radius of
#-- the fully normalized associated Legendre's function for constant order m
def clenshaw_ds_m_dr(t, q, m, clm1, slm1, lmax):
    #-- allocate for output matrix
    N = len(t)
    s_m = np.zeros((N,2),dtype=np.longdouble)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.longdouble)
    slm = scalef*slm1.astype(np.longdouble)
    #-- convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if (m == lmax):
        s_m[:,0] = -(lm + 1.0)*(clm[lmax,lmax])
        s_m[:,1] = -(lm + 1.0)*(slm[lmax,lmax])
    elif (m == (lmax-1)):
        ds_mm_dr_c_pre_1 = -(lm + 1.0) * clm[lmax,lmax-1]
        ds_mm_dr_s_pre_1 = -(lm + 1.0) * slm[lmax,lmax-1]
        a_lm = t*q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        s_m[:,0] = a_lm*ds_mm_dr_c_pre_1*clm[lmax,lmax-1]- lm*clm[lmax-1,lmax-1]
        s_m[:,1] = a_lm*ds_mm_dr_s_pre_1*slm[lmax,lmax-1]- lm*slm[lmax-1,lmax-1]
    elif ((m <= (lmax-2)) and (m >= 1)):
        ds_mm_dr_c_pre_2 = -(lm + 1.0) * clm[lmax,m]
        ds_mm_dr_s_pre_2 = -(lm + 1.0) * slm[lmax,m]
        a_lm = t*q*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        ds_mm_dr_c_pre_1 = a_lm*ds_mm_dr_c_pre_2 - lm*clm[lmax-1,m]
        ds_mm_dr_s_pre_1 = a_lm*ds_mm_dr_s_pre_2 - lm*slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*t*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))*q**2
            ds_mm_dr_c = a_lm*ds_mm_dr_c_pre_1-b_lm*ds_mm_dr_c_pre_2-(ll+1.)*clm[l,m]
            ds_mm_dr_s = a_lm*ds_mm_dr_s_pre_1-b_lm*ds_mm_dr_s_pre_2-(ll+1.)*slm[l,m]
            ds_mm_dr_c_pre_2 = np.copy(ds_mm_dr_c_pre_1)
            ds_mm_dr_s_pre_2 = np.copy(ds_mm_dr_s_pre_1)
            ds_mm_dr_c_pre_1 = np.copy(ds_mm_dr_c)
            ds_mm_dr_s_pre_1 = np.copy(ds_mm_dr_s)
        s_m[:,0] = np.copy(ds_mm_dr_c)
        s_m[:,1] = np.copy(ds_mm_dr_s)
    elif (m == 0):
        ds_mm_dr_c_pre_2 = -(lm + 1.0) * clm[lmax,0]
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))*t*q
        ds_mm_dr_c_pre_1 = a_lm * ds_mm_dr_c_pre_2 - lm * clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1)*(ll+1)))*t*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+1.)*(ll+1.))/((ll+2)*(ll+2)*(2.*ll+1.)))*q**2
            ds_mm_dr_c=a_lm*ds_mm_dr_c_pre_1-b_lm*ds_mm_dr_c_pre_2-(ll+1.)*clm[l,0]
            ds_mm_dr_c_pre_2 = np.copy(ds_mm_dr_c_pre_1)
            ds_mm_dr_c_pre_1 = np.copy(ds_mm_dr_c)
        s_m[:,0] = np.copy(ds_mm_dr_c)
    #-- return s_m rescaled with scalef
    return s_m/scalef
