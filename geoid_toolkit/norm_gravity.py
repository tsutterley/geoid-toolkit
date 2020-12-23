#!/usr/bin/env python
u"""
norm_gravity.py
Written by Tyler Sutterley (12/2020)
Calculates the normal gravity of an ellipsoid at a given latitude and height
    and calculates the derivative with respect to height

CALLING SEQUENCE:
    gamma_h, dgamma_dh = norm_gravity(latitude, height, 'WGS84')

INPUT:
    latitude: latitude in degrees
    height: height above reference ellipsoid in meters
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
    gamma_h: normal gravity for ellipsoid at height
    dgamma_dh: derivative of normal gravity with respect to height

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid

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
    Updated 12/2020: updated comments and reorganized functions
    Updated 11/2020: added function docstrings
    Updated 07/2017: added header text. higher order expansion of normal gravity
    Written 08/2013
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

def norm_gravity(latitude, height, refell):
    """
    Calculates the normal gravity of an ellipsoid at latitudes and heights
        and calculates the derivative with respect to height

    Arguments
    ---------
    latitude: latitude in degrees
    height: height above reference ellipsoid in meters
    refell: reference ellipsoid name

    Returns
    -------
    gamma_h: normal gravity for ellipsoid at height
    dgamma_dh: derivative of normal gravity with respect to height
    """

    #-- convert latitude from degrees to radians
    phi = np.pi*latitude/180.0

    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = ellip['a']
    b = ellip['b']
    #-- eccentricity
    ecc2 = ellip['ecc2']
    GM = ellip['GM']
    #-- m parameter [omega^2*a^2*b/(GM)]
    m = ellip['mp']
    #-- flattening components
    f = ellip['f']
    f_2 = -f + (5.0/2.0)*m + (1.0/2.0)*f**2.0 - \
        (26.0/7.0)*f*m + (15.0/4.0)*m**2.0
    f_4 = -(1.0/2.0)*f**2.0 + (5.0/2.0)*f*m

    #-- Normal gravity at the equator.
    #-- p. 79, Eqn.(2-186)
    gamma_a = (GM/(a * b)) * (1.0 - (3.0 / 2.0)*m - (3.0 / 14.0)*ecc2**2.0*m)
    #-- Normal gravity
    #-- p. 80, Eqn.(2-199)
    gamma_0 = gamma_a * (1.0 + f_2*np.sin(phi)**2.0 + f_4*np.sin(phi)**4.0)
    #-- Normal gravity at height
    #-- p. 82, Eqn.(2-215)
    p_1 = (1.0 + f + m - 2.0*f*np.sin(phi)**2.0)
    gamma_h = gamma_0 * (1.0 - (2.0/a)*p_1*height + (3.0/(a**2.0))*height**2.0)
    #-- approximate derivative of normal gravity with respect to height
    dgamma_dh = ((-2.0 * gamma_0) / a) * p_1
    #-- return the normal gravity and the derivative
    return (gamma_h,dgamma_dh)
