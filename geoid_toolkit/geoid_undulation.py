#!/usr/bin/env python
u"""
geoid_undulation.py
Written by Tyler Sutterley (04/2022)
Calculates the geoidal undulation at a given latitude and longitude using an
    iterative approach described in Barthelmes (2009) and Moazezi (2012)

CALLING SEQUENCE:
    geoid = geoid_undulation(lat, lon, 'WGS84', clm, slm, lmax, R, GM)

INPUT:
    lat: latitude in degrees
    lon: latitude in degrees
    refell: reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967
        GRS80 = Geodetic Reference System 1980
        WGS72 = World Geodetic System 1972
        WGS84 = World Geodetic System 1984
        ATS77 = Quasi-earth centred ellipsoid for ATS77
        NAD27 = North American Datum 1927
        NAD83 = North American Datum 1983
        INTER = International
        KRASS = Krassovsky (USSR)
        MAIRY = Modified Airy (Ireland 1965/1975)
        TOPEX = TOPEX/POSEIDON ellipsoid
        EGM96 = EGM 1996 gravity model
    clm: cosine spherical harmonics for a gravity model
    slm: sine spherical harmonics for a gravity model
    lmax: maximum spherical harmonic degree
    R: average radius used in gravity model
    GM: geocentric gravitational constant used in gravity model

OPTIONS:
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)
    EPS: level of precision for calculating geoid height

OUTPUT:
    N: geoidal undulation for a given ellipsoid in meters

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    real_potential.py: real potential at lat, lon and height for gravity model
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
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.real_potential import real_potential
from geoid_toolkit.norm_potential import norm_potential
from geoid_toolkit.norm_gravity import norm_gravity

def geoid_undulation(lat, lon, refell, clm, slm, lmax, R, GM, GAUSS=0, EPS=1e-8):
    """
    Calculates the geoidal undulation at a given latitude and longitude using
    an iterative approach described in [Barthelmes2013]_ and [Moazezi2012]_

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
    EPS: float, default 1e-8
        level of precision for calculating geoid height

    Returns
    -------
    N: float
        geoidal undulation for a given ellipsoid in meters

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

    #-- calculate the real and normal potentials for the first iteration
    W,dWdr = real_potential(lat,lon,0.0,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
    U,dUdr,dUdt = norm_potential(lat, lon, 0.0, refell, lmax)
    #-- normal gravity at latitude
    gamma_h,dgamma_dh = norm_gravity(lat, 0.0, refell)
    #-- geoid height for first iteration
    N_1 = (W - U) / gamma_h
    #-- set geoid height to the first iteration and set RMS as infinite
    N = np.copy(N_1)
    RMS = np.inf
    while (RMS > EPS):
        #-- calculate the real potentials for the iteration
        W,dWdr=real_potential(lat,lon,N_1,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
        #-- add geoid height for iteration
        N_1 += (W - U) / gamma_h
        #-- calculate RMS between iterations
        RMS = np.sqrt(np.sum((N - N_1)**2)/len(lat))
        #-- set N to the previous iteration
        N = np.copy(N_1)
    #-- return the geoid height
    return N
