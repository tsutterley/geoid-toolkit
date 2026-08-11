#!/usr/bin/env python
"""
compute.py
Written by Tyler Sutterley (04/2022)
Utilities for computing functionals from a gravity model

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    read_ICGEM_harmonics.py: reads the coefficients for a given gravity model file

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: add more keyword options to match read ICGEM options
    Updated 09/2021: define int/float precision to prevent deprecation warning
    Updated 11/2020: added function docstrings
    Updated 07/2019: split read and wrapper function into separate files
    Written 07/2017
"""

from __future__ import annotations

import pathlib
import numpy as np
import geoid_toolkit.datum
from geoid_toolkit.spatial import to_cartesian
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from geoid_toolkit.read_topography_harmonics import read_topography_harmonics


__all__ = [
    'geoid_height',
    'geoid_undulation',
    'corrected_geoid_undulation',
    'gravity_anomaly',
    'gravity_disturbance',
    'height_anomaly',
    'real_potential',
    'topographic_potential',
]


# PURPOSE: calculate geoid heights at a set of latitudes and longitudes
def geoid_height(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    gravity_model_file: str | pathlib.Path,
    **kwargs,
):
    """
    Compute geoid undulations from a gravity model

    Parameters
    ----------
    lon: float
        longitudinal points to calculate geoid height
    lat: float
        latitudinal points to calculate geoid height
    gravity_model_file: str
        full path to static gravity model file
    LMAX: int or NoneType, default None
        maximum spherical harmonic degree
    ELLIPSOID: str, default 'WGS84'
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
    TIDE: str, default 'tide_free'
        Permanent tide system of output geoid

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    GAUSS: int, default 0
        Gaussian Smoothing Radius in km
    EPS: float, default 1e-8
        level of precision for calculating geoid height
    ZIP: bool, default False
        Gravity field file is compressed in an archive file

    Returns
    -------
    N: float
        geoidal undulation for a given ellipsoid in meters
    """
    # set default keyword arguments
    kwargs.setdefault('LMAX', None)
    kwargs.setdefault('ELLIPSOID', 'WGS84')
    kwargs.setdefault('TIDE', 'tide_free')
    kwargs.setdefault('GAUSS', 0)
    kwargs.setdefault('EPS', 1e-8)
    kwargs.setdefault('ZIP', False)
    # read gravity model Ylms and change tide if specified
    Ylms = read_ICGEM_harmonics(gravity_model_file, **kwargs)
    R = np.float64(Ylms['radius'])
    GM = np.float64(Ylms['earth_gravity_constant'])
    LMAX = np.int64(Ylms['max_degree'])
    # check if reading a topography file
    # if so calculate the topographically corrected geoid undulation
    TOPOGRAPHY = kwargs.get('TOPOGRAPHY', None)
    if TOPOGRAPHY is not None:
        # read topography model Ylms
        topoYlms = read_topography_harmonics(TOPOGRAPHY, **kwargs)
        # calculate corrected geoid at coordinates
        N = corrected_geoid_undulation(
            lat,
            lon,
            kwargs['ELLIPSOID'],
            Ylms['clm'],
            Ylms['slm'],
            topoYlms['clm'],
            topoYlms['slm'],
            LMAX,
            R,
            GM,
            topoYlms['density'],
            GAUSS=kwargs['GAUSS'],
            EPS=kwargs['EPS'],
        )
    else:
        # calculate geoid at coordinates
        N = geoid_undulation(
            lat,
            lon,
            kwargs['ELLIPSOID'],
            Ylms['clm'],
            Ylms['slm'],
            LMAX,
            R,
            GM,
            GAUSS=kwargs['GAUSS'],
            EPS=kwargs['EPS'],
        )
    # return the geoid undulation
    return N


def geoid_undulation(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    GAUSS: float = 0,
    EPS: float = 1e-8,
):
    """
    Calculates the geoidal undulation using the iterative approach described in
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb`

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
    """

    # calculate the real and normal potentials for the first iteration
    W, _ = real_potential(
        lat, lon, 0.0, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
    )
    U, _, _ = geoid_toolkit.datum.norm_potential(lat, lon, 0.0, refell, lmax)
    # normal gravity at latitude
    gamma_h, _ = geoid_toolkit.datum.norm_gravity(lat, 0.0, refell)
    # geoid height for first iteration
    N_1 = (W - U) / gamma_h
    # set geoid height to the first iteration and set RMS as infinite
    N = np.copy(N_1)
    RMS = np.inf
    while RMS > EPS:
        # calculate the real potentials for the iteration
        W, _ = real_potential(
            lat, lon, N_1, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        # add geoid height for iteration
        N_1 += (W - U) / gamma_h
        # calculate RMS between iterations
        RMS = np.sqrt(np.sum((N - N_1) ** 2) / len(lat))
        # set N to the previous iteration
        N = np.copy(N_1)
    # return the geoid height
    return N


def corrected_geoid_undulation(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    tclm: float | np.ndarray,
    tslm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    density: float,
    GAUSS: float = 0,
    EPS: float = 1e-8,
):
    """
    Calculates the topographically corrected geoidal undulation
    using the iterative approach described in
    :cite:t:`Barthelmes:2013fy,Moazezi:2012fb`

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
    tclm: float
        cosine spherical harmonics for a topographic model
    tslm: float
        sine spherical harmonics for a topographic model
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
    """

    # calculate the real and normal potentials for the first iteration
    W, _ = real_potential(
        lat, lon, 0.0, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
    )
    U, _, _ = geoid_toolkit.datum.norm_potential(lat, lon, 0.0, refell, lmax)
    # topographic potential correction
    T = topographic_potential(
        lat, lon, refell, tclm, tslm, lmax, R, density, GAUSS=GAUSS
    )
    # normal gravity at latitude
    gamma_h, dgamma_dh = geoid_toolkit.datum.norm_gravity(lat, 0.0, refell)
    # geoid height for first iteration
    N_1 = (W - U - T) / gamma_h
    # set geoid height to the first iteration and set RMS as infinite
    N = np.copy(N_1)
    RMS = np.inf
    while RMS > EPS:
        # calculate the real potentials for the iteration
        W, _ = real_potential(
            lat, lon, N_1, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        # add geoid height for iteration
        N_1 += (W - U - T) / gamma_h
        # calculate RMS between iterations
        RMS = np.sqrt(np.sum((N - N_1) ** 2) / len(lat))
        # set N to the previous iteration
        N = np.copy(N_1)
    # return the geoid height
    return N


def gravity_anomaly(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    h: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    METHOD: str = 'first',
    GAUSS: int | float = 0,
):
    """
    Calculates the gravity anomaly for a given method following
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb,
    Molodensky:1958jv`

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
            - ``'molodensky'``: Molodensky method :cite:p:`Molodensky:1958jv`
    GAUSS: float, default 0
        Gaussian Smoothing Radius in km

    Returns
    -------
    ddelta_g: float
        gravity anomaly for a given ellipsoid in meters
    """
    # compute the gravity disturbance and the normal gravity
    delta_g_h = gravity_disturbance(
        lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
    )
    gamma_h, dgamma_dh = geoid_toolkit.datum.norm_gravity(lat, h, refell)
    # compute the gravity anomaly for a given method
    if METHOD.lower() == 'first':
        N = geoid_undulation(
            lat, lon, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        ddelta_g = delta_g_h + N * dgamma_dh
    elif METHOD.lower() == 'second':
        N = geoid_undulation(
            lat, lon, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        gamma_0, dgamma_d0 = geoid_toolkit.datum.norm_gravity(lat, 0, refell)
        ddelta_g = delta_g_h - (h - N) * dgamma_dh - gamma_0
    elif METHOD.lower() == 'molodensky':
        zeta = height_anomaly(
            lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        ddelta_g = delta_g_h + zeta * dgamma_dh
    return ddelta_g


def gravity_disturbance(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    h: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    GAUSS: int | float = 0,
):
    """
    Calculates the gravity disturbance following
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb,
    Molodensky:1958jv`

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
    delta_g_h: float
        gravity disturbance for a given ellipsoid in meters
    """
    # calculate the real and normal potentials at a given height
    W, dWdr = real_potential(
        lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
    )
    U, dUdr, dUdt = geoid_toolkit.datum.norm_potential(
        lat, lon, h, refell, lmax
    )
    delta_g_h = -(dWdr - dUdr)
    return delta_g_h


def height_anomaly(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    h: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    GAUSS: int | float = 0,
    EPS: float = 1e-8,
):
    """
    Calculates the height anomaly using the iterative approach described in
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb,
    Molodensky:1958jv`

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
    EPS: float, default 1e-8
        level of precision for calculating height anomaly

    Returns
    -------
    zeta: float
        height anomaly for a given ellipsoid in meters
    """
    # calculate the real and normal potentials for the first iteration
    W, dWdr = real_potential(
        lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
    )
    U, dUdr, dUdt = geoid_toolkit.datum.norm_potential(
        lat, lon, h, refell, lmax
    )
    # normal gravity at latitude
    gamma_h, dgamma_dh = geoid_toolkit.datum.norm_gravity(lat, h, refell)
    # height anomaly for first iteration
    zeta_1 = (W - U) / gamma_h
    # set zeta to the first iteration and set RMS as infinite
    zeta = np.copy(zeta_1)
    RMS = np.inf
    while RMS > EPS:
        # calculate the real and normal potentials for the iteration
        W, dWdr = real_potential(
            lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=GAUSS
        )
        U, dUdr, dUdt = geoid_toolkit.datum.norm_potential(
            lat, lon, h - zeta_1, refell, lmax
        )
        # normal gravity at latitude
        gamma_h, dgamma_dh = geoid_toolkit.datum.norm_gravity(
            lat, h - zeta_1, refell
        )
        # add height anomaly for iteration
        zeta_1 += (W - U) / gamma_h
        # calculate RMS between iterations
        RMS = np.sqrt(np.sum((zeta - zeta_1) ** 2) / len(lat))
        # set zeta to the previous iteration
        zeta = np.copy(zeta_1)
    # return the height anomaly
    return zeta


def real_potential(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    h: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    GM: float,
    GAUSS: int | float = 0,
):
    """
    Calculates the real potential using gravity model coefficients following
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb,
    Molodensky:1958jv`

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
    """

    # get ellipsoid parameters for refell
    ellip = geoid_toolkit.datum.ref_ellipsoid(refell)
    # convert coordinates to cartesian
    X, Y, Z = to_cartesian(
        lon,
        lat,
        h,
        a_axis=ellip['a'],
        flat=ellip['f'],
    )
    # height of the observation point above the ellipsoid
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # longitude and colatitude in radians
    phi = np.radians(lon)
    theta = np.pi / 2.0 - np.arctan(Z / np.hypot(X, Y))
    # number of observations
    nlat = len(lat)
    # cos and sin of colatitude
    t = np.cos(theta)
    u = np.sin(theta)
    # radius ratio
    q = (R / rr).astype(np.longdouble)

    # convert harmonics to complex form
    Ylm1 = clm - 1j * slm

    # smooth the global gravity field with a Gaussian function
    if GAUSS != 0:
        wt = 2.0 * np.pi * geoid_toolkit.math.gauss_weights(GAUSS, lmax)
        Ylm1 = np.einsum('l...,lm...->lm...', wt, Ylm1)

    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    mm = np.arange(lmax + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->pm...', mm, phi))
    # initate summations
    s_m = 0.0
    ds_m_dr = 0.0
    # iterate to calculate complete summation
    for m in range(lmax, 0, -1):
        # calculate clenshaw conditioned arrays
        cs_m = geoid_toolkit.math.clenshaw_s_m(t, q, m, Ylm1, lmax)
        dcs_m_dr = geoid_toolkit.math.clenshaw_ds_m_dr(t, q, m, Ylm1, lmax)
        # update summations and discard imaginary components
        a_m = np.sqrt((2.0 * m + 3.0) / (2.0 * m + 2.0))
        s_m *= a_m * u * q
        s_m += (cs_m * m_phi[:, m]).real
        ds_m_dr *= a_m * u * q
        ds_m_dr += (dcs_m_dr * m_phi[:, m]).real
        del cs_m, dcs_m_dr
    # calculate clenshaw conditioned arrays for order 0
    cs_m = geoid_toolkit.math.clenshaw_s_m(t, q, 0, Ylm1, lmax)
    dcs_m_dr = geoid_toolkit.math.clenshaw_ds_m_dr(t, q, 0, Ylm1, lmax)
    # add the final terms and discard imaginary components
    s_m *= np.sqrt(3.0) * u * q
    s_m += cs_m.real
    ds_m_dr *= np.sqrt(3.0) * u * q
    ds_m_dr += dcs_m_dr.real
    # compute the real potential and derivatives
    W = (GM / rr) * s_m
    dW_dr = (GM / (rr**2.0)) * ds_m_dr
    # return potentials
    return (W, dW_dr)


def topographic_potential(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    refell: str,
    clm: float | np.ndarray,
    slm: float | np.ndarray,
    lmax: int,
    R: float,
    density: float,
    GAUSS: int | float = 0,
):
    """
    Calculates the potential coefficients from a topographic model following
    :cite:t:`Barthelmes:2013fy`

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
    """
    # get ellipsoid parameters for refell
    ellip = geoid_toolkit.datum.ref_ellipsoid(refell)
    # universal gravitational constant
    G = 6.67408e-11

    # convert coordinates to cartesian
    X, Y, Z = to_cartesian(
        lon,
        lat,
        0.0,
        a_axis=ellip['a'],
        flat=ellip['f'],
    )
    # longitude and colatitude in radians
    phi = np.radians(lon)
    theta = np.pi / 2.0 - np.arctan(Z / np.hypot(X, Y))
    # number of observations
    nlat = len(lat)
    # cos and sin of colatitude
    t = np.cos(theta)
    u = np.sin(theta)
    q = 1.0

    # convert harmonics to complex form
    Ylm1 = clm - 1j * slm

    # smooth the global gravity field with a Gaussian function
    if GAUSS != 0:
        wt = 2.0 * np.pi * geoid_toolkit.math.gauss_weights(GAUSS, lmax)
        Ylm1 = np.einsum('l...,lm...->lm...', wt, Ylm1)

    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    mm = np.arange(lmax + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->pm...', mm, phi))

    # initiate summation
    s_m = 0.0
    # iterate to calculate complete summation
    for m in range(lmax, 0, -1):
        # calculate clenshaw conditioned arrays
        cs_m = geoid_toolkit.math.clenshaw_s_m(t, q, m, Ylm1, lmax)
        # update summations and discard imaginary components
        a_m = np.sqrt((2.0 * m + 3.0) / (2.0 * m + 2.0))
        s_m *= a_m * u * q
        s_m += (cs_m * m_phi[:, m]).real
        del cs_m
    # calculate clenshaw conditioned arrays for order 0
    cs_m = geoid_toolkit.math.clenshaw_s_m(t, q, 0, Ylm1, lmax)
    # add the final terms and discard imaginary components
    s_m *= np.sqrt(3.0) * u * q
    s_m += cs_m.real
    # compute the topographic potential
    T = 2.0 * np.pi * G * density * (R * s_m) ** 2
    # return the topographic potential
    return T
