#!/usr/bin/env python
u"""
gravity_anomaly.py
Written by Tyler Sutterley (04/2022)
Calculates the gravity anomaly at a given latitude and longitude using
    different methods

CALLING SEQUENCE:
    ddelta_g = gravity_anomaly(lat,lon,'WGS84',clm,slm,lmax,R,GM,METHOD='first')

INPUT:
    latitude: latitude in degrees
    longitude: latitude in degrees
    h: ellipsoidal height
    refell: reference ellipsoid specified in ref_ellipsoid.py
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

OUTPUT:
    gravity anomaly for a given ellipsoid in meters

OPTIONS:
    METHOD: method for calculating gravity anomalies
        first: classic first approximation method
        second: classic second approximation method
        molodensky: Molodensky method (1958)
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    geoid_undulation.py: geoidal undulation at a given latitude and longitude
    height_anomaly.py: height anomaly at a given latitude and longitude
    gravity_disturbance.py: gravity disturbance at a latitude and longitude
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

REFERENCE:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4
    Barthelmes, "Definition of Functionals of the Geopotential and Their
        Calculation from Spherical Harmonic Models", STR09/02 (2009)
        http://icgem.gfz-potsdam.de/str-0902-revised.pdf
    Molodensky, "New methods of studying the figure of the Earth",
        Bulletin geodesique, 50, 17-21, (1958)
        https://doi.org/10.1007/BF02537957
    Moazezi and Zomorrodian, "GGMCalc a software for calculation of the geoid
        undulation and the height anomaly using the iteration method, and
        classical gravity anomaly", Earth Science Informatics (2012)
        https://doi.org/10.1007/s12145-012-0102-2

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 01/2021: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
    Written 07/2017
"""
from geoid_toolkit.geoid_undulation import geoid_undulation
from geoid_toolkit.height_anomaly import height_anomaly
from geoid_toolkit.gravity_disturbance import gravity_disturbance
from geoid_toolkit.norm_gravity import norm_gravity

def gravity_anomaly(lat,lon,h,refell,clm,slm,lmax,R,GM,METHOD='first',GAUSS=0):
    """
    Calculates the gravity anomaly at a given latitude and longitude
    for a given method following [Barthelmes2013]_ and [Moazezi2012]_

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
    METHOD: str
        Method for calculating gravity anomalies

            - ``'first'``: classic first approximation method
            - ``'second'``: classic second approximation method
            - ``'molodensky'``: Molodensky method [Molodensky1958]_
    GAUSS: float, default 0
        Gaussian Smoothing Radius in km

    Returns
    -------
    ddelta_g: float
        gravity anomaly for a given ellipsoid in meters

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
    .. [Molodensky1958] M. S. Molodensky, "New methods of studying the figure
        of the Earth", *Bulletin G|eacute|od|eacute|sique*, 50, 17--21, (1958).
        https://doi.org/10.1007/BF02537957

    .. |eacute|    unicode:: U+00E9 .. LATIN SMALL LETTER E WITH ACUTE
    """
    #-- compute the gravity disturbance and the normal gravity
    delta_g_h = gravity_disturbance(lat,lon,h,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
    gamma_h,dgamma_dh = norm_gravity(lat, h, refell)
    #-- compute the gravity anomaly for a given method
    if (METHOD.lower() == 'first'):
        N = geoid_undulation(lat,lon,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
        ddelta_g = delta_g_h + N * dgamma_dh
    elif (METHOD.lower() == 'second'):
        N = geoid_undulation(lat,lon,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
        gamma_0,dgamma_d0 = norm_gravity(lat, 0, refell)
        ddelta_g = delta_g_h - (h-N) * dgamma_dh - gamma_0
    elif (METHOD.lower() == 'molodensky'):
        zeta = height_anomaly(lat,lon,h,refell,clm,slm,lmax,R,GM,GAUSS=GAUSS)
        ddelta_g = delta_g_h + zeta * dgamma_dh
    return ddelta_g
