#!/usr/bin/env python
"""
topographic_potential.py
Written by Tyler Sutterley (07/2026)
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
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use complex form of spherical harmonics for summations
        use np.radians to convert from degrees to radians
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2017
"""

import numpy as np
from geoid_toolkit.spatial import to_cartesian
from geoid_toolkit.gauss_weights import gauss_weights
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid


def topographic_potential(
    lat, lon, refell, clm, slm, lmax, R, density, GAUSS=0
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
    ellip = ref_ellipsoid(refell)
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

    # convert harmonics to complex form
    Ylm1 = clm - 1j * slm

    # smooth the global gravity field with a Gaussian function
    if GAUSS != 0:
        wt = 2.0 * np.pi * gauss_weights(GAUSS, lmax)
        Ylm1 = np.einsum('l...,lm...->lm...', wt, Ylm1)

    # calculate clenshaw summations
    cs_m = np.zeros((nlat, lmax + 1), dtype=np.clongdouble)
    for m in range(lmax, -1, -1):
        cs_m[:, m] = _clenshaw_s_m(t, m, Ylm1, lmax)

    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    mm = np.arange(lmax + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->pm...', mm, phi))

    # calculate summation and drop imaginary component
    s_m = (cs_m[:, lmax] * m_phi[:, lmax]).real
    # iterate to calculate complete summation
    for m in range(lmax - 1, 0, -1):
        # update summations and discard imaginary components
        a_m = np.sqrt((2.0 * m + 3.0) / (2.0 * m + 2.0))
        s_m = a_m * u * s_m + (cs_m[:, m] * m_phi[:, m]).real
    # add the final terms
    s_m = np.sqrt(3.0) * u * s_m + cs_m[:, 0].real
    # compute the topographic potential
    T = 2.0 * np.pi * G * density * (R * s_m) ** 2
    # return the topographic potential
    return T


# PURPOSE: compute Clenshaw summation of the fully normalized associated
# Legendre's function for constant order m
def _clenshaw_s_m(t, m, Ylm1, lmax, SCALE=1e-280):
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
        a_lm = t * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        cs_m[:] = a_lm * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = t * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = t * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
            )
            b_lm = np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:] = np.copy(s_mm_l)
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = t * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = t * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0)) / ((ll + 1) * (ll + 1))
            )
            b_lm = np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:] = np.copy(s_mm_l)
    # return rescaled cs_m
    return cs_m / SCALE
