#!/usr/bin/env python
u"""
height_anomaly.py
Written by Tyler Sutterley (01/2021)
Calculates the height anomaly at a given latitude and longitude using an
    iterative approach described in Barthelmes (2009) and Moazezi (2012)

CALLING SEQUENCE:
    ha = height_anomaly(lat, lon, 'WGS84', clm, slm, lmax, R, GM)

INPUT:
    latitude: latitude in degrees
    longitude: latitude in degrees
    h: ellipsoidal height
    refell: reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967
        GRS80 = Geodetic Reference System 1980
        WGS72 = World Geodetic System 1972
        WGS84 = World Geodetic System 1984
        ATS77 = Quasi-earth centred ellipsoid for ATS77
        NAD27 = North American Datum 1927 (=CLK66)
        NAD83 = North American Datum 1983 (=GRS80)
        INTER = International
        KRASS = Krassovsky (USSR)
        MAIRY = Modified Airy (Ireland 1965/1975)
        TOPEX = TOPEX/POSEIDON ellipsoid
        EGM96 = EGM 1996 gravity model
    clm: cosine spherical harmonics for a gravity model
    slm: sine spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model

OUTPUT:
     height anomaly for a given ellipsoid in meters

OPTIONS:
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)
    EPS: level of precision for calculating height anomaly

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
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

UPDATE HISTORY:
    Updated 01/2021: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.real_potential import real_potential
from geoid_toolkit.norm_potential import norm_potential
from geoid_toolkit.norm_gravity import norm_gravity

def height_anomaly(lat,lon,h,refell,clm,slm,lmax,R,GM,GAUSS=0,EPS=1e-8):
    """
    Calculates the height anomaly at a given latitude and longitude using an
        iterative approach

    Arguments
    ---------
    lat: latitude in degrees
    lon: latitude in degrees
    h: ellipsoidal height
    refell: reference ellipsoid name
    clm: cosine spherical harmonics for a gravity model
    slm: sine spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model

    Keyword arguments
    -----------------
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)
    EPS: level of precision for calculating height anomaly

    Returns
    -------
    zeta: height anomaly for a given ellipsoid in meters
    """
    #-- calculate the real and normal potentials for the first iteration
    W,dWdr = real_potential(lat,lon,h,refell,GM,R,clm,slm,lmax,GAUSS=GAUSS)
    U,dUdr,dUdt = norm_potential(lat, lon, h, refell, lmax)
    #-- normal gravity at latitude
    gamma_h,dgamma_dh = norm_gravity(lat, h, refell)
    #-- height anomaly for first iteration
    zeta_1 = (W - U) / gamma_h
    #-- set zeta to the first iteration and set RMS as infinite
    zeta = np.copy(zeta_1)
    RMS = np.inf
    while (RMS > EPS):
        #-- calculate the real and normal potentials for the iteration
        W,dWdr = real_potential(lat,lon,h,refell,GM,R,clm,slm,lmax,GAUSS=GAUSS)
        U,dUdr,dUdt = norm_potential(lat,lon,h-zeta_1,refell,lmax)
        #-- normal gravity at latitude
        gamma_h,dgamma_dh = norm_gravity(lat,h-zeta_1,refell)
        #-- add height anomaly for iteration
        zeta_1 += (W - U) / gamma_h
        #-- calculate RMS between iterations
        RMS = np.sqrt(np.sum((zeta - zeta_1)**2)/len(lat))
        #-- set zeta to the previous iteration
        zeta = np.copy(zeta_1)
    #-- return the height anomaly
    return zeta
