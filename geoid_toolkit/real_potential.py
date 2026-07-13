#!/usr/bin/env python
"""
real_potential.py
Written by Tyler Sutterley (07/2026)
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
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use complex form of spherical harmonics for summations
        use np.radians to convert from degrees to radians
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
        changed dtypes to long double for high degree and order models
    Written 07/2017
"""

import numpy as np
from geoid_toolkit.spatial import to_cartesian
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
from geoid_toolkit.gauss_weights import gauss_weights


def real_potential(lat, lon, h, refell, clm, slm, lmax, R, GM, GAUSS=0):
    """
    Calculates the real potential using gravity model coefficients following
    :cite:t:`Barthelmes:2013fy,HofmannWellenhof:2006hy,Moazezi:2012fb,Molodensky:1958jv`

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
    """

    # get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    # convert from geodetic latitude to geocentric latitude
    phi = np.radians(lon)
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
    # colatitude in radians
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
        wt = 2.0 * np.pi * gauss_weights(GAUSS, lmax)
        Ylm1 = np.einsum('l...,lm...->lm...', wt, Ylm1)

    # calculate clenshaw summations
    cs_m = np.zeros((nlat, lmax + 1), dtype=np.clongdouble)
    dcs_m_dr = np.zeros((nlat, lmax + 1), dtype=np.clongdouble)
    for m in range(lmax, -1, -1):
        cs_m[:, m] = _clenshaw_s_m(t, q, m, Ylm1, lmax)
        dcs_m_dr[:, m] = _clenshaw_ds_m_dr(t, q, m, Ylm1, lmax)

    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    m_phi = np.exp(1j * np.einsum("m...,p...->pm...", m, phi))
    # calculate summation and drop imaginary component
    s_m = (cs_m[:, lmax] * m_phi[:, lmax]).real
    ds_m_dr = (dcs_m_dr[:, lmax] * m_phi[:, lmax]).real

    # iterate to calculate complete summation
    for m in range(lmax - 1, 0, -1):
        # update summations and discard imaginary components
        a_m = np.sqrt((2.0 * m + 3.0) / (2.0 * m + 2.0))
        s_m = a_m * u * q * s_m + (cs_m[:, m] * m_phi[:, m]).real
        ds_m_dr = a_m * u * q * ds_m_dr + (dcs_m_dr[:, m] * m_phi[:, m]).real

    # add the final terms
    s_m = np.sqrt(3.0) * u * q * s_m + cs_m[:, 0].real
    ds_m_dr = np.sqrt(3.0) * u * q * ds_m_dr + dcs_m_dr[:, 0].real
    # compute the real potential and derivatives
    W = (GM / rr) * s_m
    dW_dr = (GM / (rr**2.0)) * ds_m_dr
    # return the real potential and derivatives
    return (W, dW_dr)


# PURPOSE: compute Clenshaw summation of the fully normalized associated
# Legendre's function for constant order m
def _clenshaw_s_m(t, q, m, Ylm1, lmax, SCALE=1e-280):
    # allocate for output matrix
    N = len(t)
    cs_m = np.zeros((N), dtype=np.clongdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.clongdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        cs_m[:] = np.copy(ylm[lmax, lmax])
    elif m == (lmax - 1):
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        cs_m[:] = a_lm * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:] = np.copy(s_mm_l)
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = (
            t * q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1) * (ll + 1))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:, 0] = np.copy(s_mm_l)
    # return rescaled cs_m
    return cs_m / SCALE


# PURPOSE: compute Clenshaw summation of derivative with respect to latitude
# of the fully normalized associated Legendre's function for constant order m
def _clenshaw_ds_m(t, u, q, m, Ylm1, lmax, SCALE=1e-280):
    # allocate for output matrix
    N = len(t)
    dcs_m = np.zeros((N, 2), dtype=np.longdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.longdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        dcs_m[:] = mm * t * u * ylm[lmax, lmax]
    elif m == (lmax - 1):
        a_lm = q * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        s_dot_mm = a_lm * ylm[lmax, lmax - 1]
        s_mm_l = a_lm * t * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
        dcs_m[:] = mm * t * u * s_mm_l - u * s_dot_mm
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = q * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        s_dot_mm_minus_2 = 0.0
        s_dot_mm_minus_1 = a_lm * s_mm_minus_2
        s_mm_minus_1 = a_lm * t * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = q * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_dot_mm = (
                a_lm * (s_dot_mm_minus_1 * t + s_mm_minus_1)
                - b_lm * s_dot_mm_minus_2
            )
            s_mm_l = a_lm * t * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
            s_dot_mm_minus_2 = np.copy(s_dot_mm_minus_1)
            s_dot_mm_minus_1 = np.copy(s_dot_mm)
        dcs_m[:] = mm * t * u * s_mm_l - u * s_dot_mm
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        s_dot_mm_minus_1 = a_lm * s_mm_minus_2
        s_mm_minus_1 = a_lm * t * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = q * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0)) / ((ll + 1) * (ll + 1))
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            s_dot_mm = (
                a_lm * (s_dot_mm_minus_1 * t + s_mm_minus_1)
                - b_lm * s_dot_mm_minus_2
            )
            s_mm_l = a_lm * t * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
            s_dot_mm_minus_2 = np.copy(s_dot_mm_minus_1)
            s_dot_mm_minus_1 = np.copy(s_dot_mm)
        dcs_m[:, 0] = -u * s_dot_mm
    # return rescaled dcs_m
    return dcs_m / SCALE


# PURPOSE: compute Clenshaw summation of derivative with respect to radius of
# the fully normalized associated Legendre's function for constant order m
def _clenshaw_ds_m_dr(t, q, m, Ylm1, lmax, SCALE=1e-280):
    # allocate for output matrix
    N = len(t)
    dcs_m_dr = np.zeros((N), dtype=np.longdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.longdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        dcs_m_dr[:] = -(lm + 1.0) * (ylm[lmax, lmax])
    elif m == (lmax - 1):
        ds_mm_dr_c_pre_1 = -(lm + 1.0) * ylm[lmax, lmax - 1]
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        dcs_m_dr[:] = (
            a_lm * ds_mm_dr_c_pre_1 * ylm[lmax, lmax - 1]
            - lm * ylm[lmax - 1, lmax - 1]
        )
    elif (m <= (lmax - 2)) and (m >= 1):
        ds_mm_dr_minus_2 = -(lm + 1.0) * ylm[lmax, m]
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        ds_mm_dr_minus_1 = a_lm * ds_mm_dr_minus_2 - lm * ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            ds_mm_dr = (
                a_lm * ds_mm_dr_minus_1
                - b_lm * ds_mm_dr_minus_2
                - (ll + 1.0) * ylm[l, m]
            )
            ds_mm_dr_minus_2 = np.copy(ds_mm_dr_minus_1)
            ds_mm_dr_minus_1 = np.copy(ds_mm_dr)
        dcs_m_dr[:] = np.copy(ds_mm_dr)
    elif m == 0:
        ds_mm_dr_minus_2 = -(lm + 1.0) * ylm[lmax, 0]
        a_lm = (
            t * q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        )
        ds_mm_dr_minus_1 = a_lm * ds_mm_dr_minus_2 - lm * ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1) * (ll + 1))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            ds_mm_dr = (
                a_lm * ds_mm_dr_minus_1
                - b_lm * ds_mm_dr_minus_2
                - (ll + 1.0) * ylm[l, 0]
            )
            ds_mm_dr_minus_2 = np.copy(ds_mm_dr_minus_1)
            ds_mm_dr_minus_1 = np.copy(ds_mm_dr)
        dcs_m_dr[:] = np.copy(ds_mm_dr)
    # return rescaled dcs_m_dr
    return dcs_m_dr / SCALE
