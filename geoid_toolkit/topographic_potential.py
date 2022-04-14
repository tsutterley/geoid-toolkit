#!/usr/bin/env python
u"""
topographic_potential.py
Written by Tyler Sutterley (04/2022)
Calculates the potential at a given latitude and height using
    coefficients from a topographic model

CALLING SEQUENCE:
    T = topographic_potential(lat, lon, clm, slm, lmax, R, density)

INPUT:
    latitude: latitude in degrees
    longitude: longitude in degrees
    clm: cosine spherical harmonics for a topographic model
    slm: sin spherical harmonics for a topographic model
    lmax: maximum spherical harmonic degree
    R: mean radius of the Earth using parameters for a gravity model
    density: density of the topography in the model

OUTPUT:
    T: potential from topography model

OPTIONS:
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

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
        http://icgem.gfz-potsdam.de/ICGEM/theory/str-0902-revised.pdf
    Moazezi and Zomorrodian, "GGMCalc a software for calculation of the geoid
        undulation and the height anomaly using the iteration method, and
        classical gravity anomaly", Earth Science Informatics (2012)
        http://dx.doi.org/10.1007/s12145-012-0102-2
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        http://dx.doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.gauss_weights import gauss_weights
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

def topographic_potential(lat, lon, refell, clm, slm, lmax, R, density, GAUSS=0):
    """
    Calculates the potential at a given latitude and height using
    coefficients from a topographic model following [Barthelmes2013]_

    Parameters
    ----------
    lat: float
        latitude in degrees
    lon: float
        longitude in degrees
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
        cosine spherical harmonics for a topographic model
    slm: float
        sine spherical harmonics for a topographic model
    lmax: int
        maximum spherical harmonic degree
    R: float
        average radius used in gravity model
    density: float
        density of the topography in the model
    GAUSS: float, default 0
        Gaussian Smoothing Radius in km

    Returns
    -------
    T: float
        potential from topography model

    References
    ----------
    .. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the
        Geopotential and Their Calculation from Spherical Harmonic Models",
        *GeoForschungsZentrum Scientific Technical Report*, STR09/02, (2013).
        `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_
    """
    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = ellip['a']
    ecc1 = ellip['ecc1']
    #-- universal gravitational constant
    G = 6.67408e-11

    #-- convert from geodetic latitude to geocentric latitude
    latitude_geodetic_rad = np.pi*lat/180.0
    longitude_rad = np.pi*lon/180.0
    #-- prime vertical radius of curvature
    N = a/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    X = N * np.cos(latitude_geodetic_rad) * np.cos(longitude_rad)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    #-- number of observations
    nlat = len(lat)
    #-- sin and cos of latitude
    latitude_geocentric_rad = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))
    t = np.sin(latitude_geocentric_rad)
    u = np.cos(latitude_geocentric_rad)

    #-- smooth the global gravity field with a Gaussian function
    if (GAUSS != 0):
        wt = 2.0*np.pi*gauss_weights(GAUSS,lmax)
        for l in range(0,lmax+1):
            clm[l,:] = clm[l,:]*wt[l]
            slm[l,:] = slm[l,:]*wt[l]

    #-- calculate clenshaw summations
    s_m_c = np.zeros((nlat,lmax*2+2))
    for m in range(lmax, -1, -1):
        s_m_c[:,2*m:2*m+2] = clenshaw_s_m(t,m,clm,slm,lmax)

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
    #-- iterate to calculate complete summation
    for m in range(lmax-1, 0, -1):
        cos_m_phi[:,m] = cos_phi_2*cos_m_phi[:,m+1] - cos_m_phi[:,m+2]
        sin_m_phi[:,m] = cos_phi_2*sin_m_phi[:,m+1] - sin_m_phi[:,m+2]
        a_m = np.sqrt((2.0*m+3.0)/(2.0*m+2.0))
        s_m = a_m*u*s_m + s_m_c[:,2*m]*cos_m_phi[:,m] + s_m_c[:,2*m+1]*sin_m_phi[:,m]

    #-- compute the topographic potential
    s_m = np.sqrt(3.0)*u*s_m + s_m_c[:,0]
    T = 2.0 * np.pi * G * density * (R * s_m)**2
    #-- return the topographic potential
    return T

#-- PURPOSE: compute Clenshaw summation of the fully normalized associated
#-- Legendre's function for constant order m
def clenshaw_s_m(t, m, clm1, slm1, lmax):
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
        a_lm = t*np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))
        s_m[:,0] = a_lm*clm[lmax,lmax-1] + clm[lmax-1,lmax-1]
        s_m[:,1] = a_lm*slm[lmax,lmax-1] + slm[lmax-1,lmax-1]
    elif ((m <= (lmax-2)) and (m >= 1)):
        s_mm_c_pre_2 = np.copy(clm[lmax,m])
        s_mm_s_pre_2 = np.copy(slm[lmax,m])
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t
        s_mm_c_pre_1 = a_lm*s_mm_c_pre_2 + clm[lmax-1,m]
        s_mm_s_pre_1 = a_lm*s_mm_s_pre_2 + slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*t
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))
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
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))*t
        s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = np.longdouble(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1)*(ll+1)))*t
            b_lm=np.sqrt(((2.*ll+5.)*(ll+1.)*(ll+1.))/((ll+2)*(ll+2)*(2.*ll+1.)))
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,0]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
        s_m[:,0] = np.copy(s_mm_c)
    #-- return s_m rescaled with scalef
    return s_m/scalef
