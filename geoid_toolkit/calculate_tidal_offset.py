#!/usr/bin/env python
u"""
calculate_tidal_offset.py
Written by Tyler Sutterley (11/2020)
Calculates the spherical harmonic offset for a tide system to change from a tide
    free state where there is no permanent direct and indirect tidal potentials

CALLING SEQUENCE:
    delta = calculate_tidal_offset(TIDE, GM, R, refell)

INPUT:
    TIDE: output tidal system (changing from tide free)
        mean_tide: restores permanent tidal potentials (direct and indirect)
        zero_tide: restores permanent direct tidal potential
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model
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

OUTPUT:
    delta: offset for changing from tide free system

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
    Updated 11/2020: added function docstrings
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

def calculate_tidal_offset(TIDE, GM, R, refell):
    """
    Calculates the spherical harmonic offset for a tide system to change
        from a tide free state where there is no permanent direct and
        indirect tidal potentials

    Arguments
    ---------
    TIDE: output tidal system
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model
    refell: reference ellipsoid name

    Returns
    -------
    deltaC20: offset for changing from tide free system
    """
    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    #-- standard gravitational acceleration
    gamma = 9.80665
    trans = (-0.198*gamma*R**3)/(np.sqrt(5.0)*GM*ellip['a']**2)
    #-- load love number for degree 2 from PREM (Han and Wahr, 1995)
    k2 = -0.30252982142510
    #-- conversion for each tidal system
    if (TIDE == 'mean_tide'):
        conv = (1.0 + k2)
    elif (TIDE == 'zero_tide'):
        conv = k2
    #-- return the C20 offset
    return conv*trans
