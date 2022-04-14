#!/usr/bin/env python
u"""
calculate_tidal_offset.py
Written by Tyler Sutterley (04/2022)
Calculates the spherical harmonic offset to change tide systems

CALLING SEQUENCE:
    delta = calculate_tidal_offset(TIDE, GM, R, refell)

INPUT:
    TIDE: output tidal system
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
    R: average radius used in gravity model
    GM: geocentric gravitational constant used in gravity model
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

OPTIONS:
    LOVE: load love number for degree 2
    REFERENCE: original tidal system
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed

OUTPUT:
    delta: offset for changing to tide system

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid

REFERENCE:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4
    Losch and Seufer, "How to Compute Geoid Undulations (Geoid Height Relative
        to a Given Reference Ellipsoid) from Spherical Harmonic Coefficients for
        Satellite Altimetry Applications" (2003)
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: can change from different tidal systems and to tide free
    Updated 11/2020: added function docstrings
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

def calculate_tidal_offset(TIDE, GM, R, refell, LOVE=0.3,
    REFERENCE='tide_free'):
    """
    Calculates the spherical harmonic offset to change permanent tide systems

    Parameters
    ----------
    TIDE: str
        Output permanent tidal system

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    R: float
        Average radius used in gravity model
    GM: float
        Geocentric gravitational constant used in gravity model
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
    LOVE: float, default 0.3
        Load love number for degree 2
    REFERENCE: str, default 'tide_free'
        Original permanent tidal system of gravity modeld

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed

    Returns
    -------
    delta: float
        Offset for changing to tide system

    References
    ----------
    .. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz,
        *Physical Geodesy*, 2nd Edition, 403 pp., (2006).
        `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
    .. [Losch2003] M. Losch and V. Seufer,
        "How to Compute Geoid Undulations (Geoid Height Relative
        to a Given Reference Ellipsoid) from Spherical Harmonic
        Coefficients for Satellite Altimetry Applications", (2003).
        `eprint ID: 11802 <http://mitgcm.org/~mlosch/geoidcookbook.pdf>`_
    """
    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    #-- standard gravitational acceleration
    gamma = 9.80665
    trans = (-0.198*gamma*R**3)/(np.sqrt(5.0)*GM*ellip['a']**2)
    #-- conversion to switch to tide free
    if (REFERENCE == 'tide_free'):
        tide_free_conv = 0.0
    elif (REFERENCE == 'mean_tide'):
        tide_free_conv = -(1.0 + LOVE)
    elif (REFERENCE == 'zero_tide'):
        tide_free_conv = -LOVE
    #-- conversion for each tidal system
    if (TIDE == 'mean_tide'):
        conv = (1.0 + LOVE) + tide_free_conv
    elif (TIDE == 'zero_tide'):
        conv = LOVE + tide_free_conv
    elif (TIDE == 'tide_free'):
        conv = 0.0 + tide_free_conv
    #-- return the C20 offset to change tide systems
    delta = conv*trans
    return delta
