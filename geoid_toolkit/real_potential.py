#!/usr/bin/env python
u"""
real_potential.py
Written by Tyler Sutterley (11/2020)
Calculates the real potential at a given latitude and height using
    coefficients from a gravity model

CALLING SEQUENCE:
    W, dW_dr, dW_dtheta = real_potential(lat, lon, h, clm, slm, lmax)

INPUT:
    latitude: latitude in degrees
    longitude: longitude in degrees
    height: height above reference ellipsoid in meters
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model
    clm: cosine spherical harmonics for a gravity model
    slm: sin spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree

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
    Updated 11/2020: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
        changed dtypes to long double for high degree and order models
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
from geoid_toolkit.gauss_weights import gauss_weights

def real_potential(latitude, longitude, h, refell, GM, R, clm, slm,
    lmax, GAUSS=0):
    """
    Calculates the real potential at a given latitude and height using
        coefficients from a gravity model

    Arguments
    ---------
    latitude: latitude in degrees
    longitude: longitude in degrees
    height: height above reference ellipsoid in meters
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model
    clm: cosine spherical harmonics for a gravity model
    slm: sin spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree

    Keyword arguments
    -----------------
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

    Returns
    -------
    W: real potential at height h
    dW_dr: derivative of real potential with respect to radius
    dW_dtheta: derivative of real potential with respect to theta
    """

    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = np.float128(ellip['a'])
    ecc1 = np.float128(ellip['ecc1'])

    #-- convert from geodetic latitude to geocentric latitude
    latitude_geodetic_rad = (np.pi*latitude/180.0).astype(np.float128)
    longitude_rad = (np.pi*longitude/180.0).astype(np.float128)
    #-- prime vertical radius of curvature
    N = a/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(longitude_rad)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    #-- height of the observation point above the ellipsoid
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    latitude_geocentric_rad = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))
    #-- number of observations
    nlat = len(latitude)
    #-- sin and cos of latitude
    t = np.sin(latitude_geocentric_rad)
    u = np.cos(latitude_geocentric_rad)
    #-- radius ratio
    q = (R / rr).astype(np.float128)

    #-- smooth the global gravity field with a Gaussian function
    if (GAUSS != 0):
        wt = 2.0*np.pi*gauss_weights(GAUSS,lmax)
        for l in range(0,lmax+1):
            clm[l,:] = clm[l,:]*wt[l]
            slm[l,:] = slm[l,:]*wt[l]

    #-- calculate clenshaw summations
    s_m_c = np.zeros((nlat,lmax*2+2),dtype=np.float128)
    ds_m_dr_c = np.zeros((nlat,lmax*2+2),dtype=np.float128)
    for m in range(lmax, -1, -1):
        s_m_c[:,2*m:2*m+2] = clenshaw_s_m(t,q,m,clm,slm,lmax)
        ds_m_dr_c[:,2*m:2*m+2] = clenshaw_ds_m_dr(t,q,m,clm,slm,lmax)

    #-- calculate cos phi
    cos_phi_2 = 2.0*np.cos(longitude_rad)
    #-- matrix of cos/sin m*phi (longitude_rad) summation
    cos_m_phi = np.zeros((nlat,lmax+2),dtype=np.float128)
    sin_m_phi = np.zeros((nlat,lmax+2),dtype=np.float128)
    #-- initialize matrix with values at lmax+1 and lmax
    cos_m_phi[:,lmax+1] = np.cos(np.float128(lmax + 1)*longitude_rad)
    sin_m_phi[:,lmax+1] = np.sin(np.float128(lmax + 1)*longitude_rad)
    cos_m_phi[:,lmax] = np.cos(np.float128(lmax)*longitude_rad)
    sin_m_phi[:,lmax] = np.sin(np.float128(lmax)*longitude_rad)
    #-- calculate summation
    s_m = s_m_c[:,2*lmax]*cos_m_phi[:,lmax] + s_m_c[:,2*lmax+1]*sin_m_phi[:,lmax]
    ds_m_dr = ds_m_dr_c[:,2*lmax]*cos_m_phi[:,lmax] + \
        ds_m_dr_c[:,2*lmax+1]*sin_m_phi[:,lmax]

    #-- iterate to calculate complete summation
    for m in range(lmax-1, 0, -1):
        cos_m_phi[:,m] = cos_phi_2*cos_m_phi[:,m+1] - cos_m_phi[:,m+2]
        sin_m_phi[:,m] = cos_phi_2*sin_m_phi[:,m+1] - sin_m_phi[:,m+2]
        a_m = np.sqrt((2.0*np.float128(m)+3.0)/(2.0*np.float128(m)+2.0))
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
    s_m = np.zeros((N,2),dtype=np.float128)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.float128)
    slm = scalef*slm1.astype(np.float128)
    #-- convert lmax and m to float
    lm = np.float128(lmax)
    mm = np.float128(m)
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
            ll = np.float128(l)
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
            ll = np.float128(l)
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
    s_m = np.zeros((N,2),dtype=np.float128)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.float128)
    slm = scalef*slm1.astype(np.float128)
    #-- convert lmax and m to float
    lm = np.float128(lmax)
    mm = np.float128(m)
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
            ll = np.float128(l)
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
            ll = np.float128(l)
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
    s_m = np.zeros((N,2),dtype=np.float128)
    #-- scaling factor to prevent overflow
    scalef = 1.0e-280
    clm = scalef*clm1.astype(np.float128)
    slm = scalef*slm1.astype(np.float128)
    #-- convert lmax and m to float
    lm = np.float128(lmax)
    mm = np.float128(m)
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
            ll = np.float128(l)
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
            ll = np.float128(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1)*(ll+1)))*t*q
            b_lm=np.sqrt(((2.*ll+5.)*(ll+1.)*(ll+1.))/((ll+2)*(ll+2)*(2.*ll+1.)))*q**2
            ds_mm_dr_c=a_lm*ds_mm_dr_c_pre_1-b_lm*ds_mm_dr_c_pre_2-(ll+1.)*clm[l,0]
            ds_mm_dr_c_pre_2 = np.copy(ds_mm_dr_c_pre_1)
            ds_mm_dr_c_pre_1 = np.copy(ds_mm_dr_c)
        s_m[:,0] = np.copy(ds_mm_dr_c)
    #-- return s_m rescaled with scalef
    return s_m/scalef
