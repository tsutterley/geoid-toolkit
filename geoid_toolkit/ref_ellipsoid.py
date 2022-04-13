#!/usr/bin/env python
u"""
ref_ellipsoid.py
Written by Tyler Sutterley (04/2022)

Computes parameters for a reference ellipsoid

CALLING SEQUENCE
    wgs84 = ref_ellipsoid('WGS84',UNITS='MKS')

INPUT:
    refell - reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967 (IAU ellipsoid)
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
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data

OPTIONS:
    UNITS: output units
        MKS: meters, kilograms, seconds
        CGS: centimeters, grams, seconds

OUTPUTS:
    a: semimajor semi-axis (m)
    b: semiminor semi-axis (m)
    f: flattening
    c: Polar radius of curvature
    rad_e: mean radius of ellipsoid having the same volume
    rad_p: Polar radius of curvature
    C20: Normalized C20 harmonic
    norm_a: Normal gravity at the equator
    norm_b: Normal gravity at the pole
    U0: Normal potential at the ellipsoid
    dk: ratio between gravity at pole versus gravity at equator
    m: m parameter (m)
    lin_ecc: Linear eccentricity
    ecc1: First eccentricity
    ecc2: Second eccentricity
    area: Area of the ellipsoid

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Hofmann-Wellenhof and Moritz (2006)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: added function docstrings
    Updated 09/2017: added parameters for Hughes (1980) ellipsoid
    Updated 07/2017: output J2 harmonic (in addition to C20)
    Updated 06/2017: added parameters for EGM1996. CGS now True/False
    Updated 06/2016: using __future__ print function, formatted output
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 04/2015: added command line option.  little updates to comments
    Updated 06/2014: changed message to sys.exit
    Updated 02/2014: minor update to if statements
"""
import numpy as np

def ref_ellipsoid(refell, UNITS='MKS'):
    """
    Computes parameters for a reference ellipsoid

    Parameters
    ----------
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
    UNITS: str, default 'MKS'
        Output units

            - ``'MKS'``: meters, kilograms, seconds
            - ``'CGS'``: centimeters, grams, seconds

    Returns
    -------
    a: float
        semimajor semi-axis (m)
    b: float
        semiminor semi-axis (m)
    f: float
        flattening
    c: float
        Polar radius of curvature
    rad_e: float
        Mean radius of ellipsoid having the same volume
    rad_p: float
        Polar radius of curvature
    C20: float
        Normalized C20 harmonic
    norm_a: float
        Normal gravity at the equator
    norm_b: float
        Normal gravity at the pole
    U0: float
        Normal potential at the ellipsoid
    dk: float
        ratio between gravity at pole versus gravity at equator
    m: float
        m parameter (m)
    lin_ecc: float
        Linear eccentricity
    ecc1: float
        First eccentricity
    ecc2: float
        Second eccentricity
    area: float
        Area of the ellipsoid

    References
    ----------
    .. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz,
        *Physical Geodesy*, 2nd Edition, 403 pp., (2006).
        `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
    """

    if refell.upper() in ('CLK66','NAD27'):
        #-- Clarke 1866
        a_axis = 6378206.4#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/294.9786982#-- flattening of the ellipsoid

    elif refell.upper() in ('GRS80','NAD83'):
        #-- Geodetic Reference System 1980
        #-- North American Datum 1983
        a_axis = 6378135.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.26#-- flattening of the ellipsoid
        GM = 3.986005e14#-- [m^3/s^2] Geocentric Gravitational Constant

    elif (refell.upper() == 'GRS67'):
        #-- Geodetic Reference System 1967
        #-- International Astronomical Union (IAU ellipsoid)
        a_axis = 6378160.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.247167427#-- flattening of the ellipsoid
        GM = 3.98603e14#-- [m^3/s^2] Geocentric Gravitational Constant
        omega = 7292115.1467e-11#-- angular velocity of the Earth [rad/s]

    elif (refell.upper() == 'WGS72'):
        #-- World Geodetic System 1972
        a_axis = 6378135.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.26#-- flattening of the ellipsoid

    elif (refell.upper() == 'WGS84'):
        #-- World Geodetic System 1984
        a_axis = 6378137.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.257223563#-- flattening of the ellipsoid

    elif (refell.upper() == 'ATS77'):
        #-- Quasi-earth centred ellipsoid for ATS77
        a_axis = 6378135.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.257#-- flattening of the ellipsoid

    elif (refell.upper() == 'KRASS'):
        #-- Krassovsky (USSR)
        a_axis = 6378245.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.3#-- flattening of the ellipsoid

    elif (refell.upper() == 'INTER'):
        #-- International
        a_axis = 6378388.0#-- [m] semimajor axis of the ellipsoid
        flat = 1/297.0#-- flattening of the ellipsoid

    elif (refell.upper() == 'MAIRY'):
        #-- Modified Airy (Ireland 1965/1975)
        a_axis = 6377340.189#-- [m] semimajor axis of the ellipsoid
        flat = 1/299.3249646#-- flattening of the ellipsoid

    elif (refell.upper() == 'TOPEX'):
        #-- TOPEX/POSEIDON ellipsoid
        a_axis = 6378136.3#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.257#-- flattening of the ellipsoid
        GM = 3.986004415e14#-- [m^3/s^2]

    elif (refell.upper() == 'EGM96'):
        #-- EGM 1996 gravity model
        a_axis = 6378136.3#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.256415099#-- flattening of the ellipsoid
        GM = 3.986004415e14#-- [m^3/s^2]

    elif (refell.upper() == 'HGH80'):
        #-- Hughes 1980 Ellipsoid used in some NSIDC data
        a_axis = 6378273.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.279411123064#-- flattening of the ellipsoid

    else:
        raise ValueError('Incorrect reference ellipsoid Name')

    if refell.upper() not in ('GRS80','GRS67','NAD83','TOPEX','EGM96'):
        #-- for ellipsoids not listing the Geocentric Gravitational Constant
        GM = 3.986004418e14#-- [m^3/s^2]

    if refell.upper() not in ('GRS67'):
        #-- for ellipsoids not listing the angular velocity of the Earth
        omega = 7292115e-11#-- [rad/s]

    #-- convert units to CGS
    if (UNITS == 'CGS'):
        a_axis *= 100.0
        GM *= 10e6

    #-- DERIVED PARAMETERS:
    #-- mean radius of the Earth having the same volume
    #-- (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1D -f)
    rad_e = a_axis*(1.0 -flat)**(1.0/3.0)

    #-- semiminor axis of the ellipsoid
    b_axis = (1.0 -flat)*a_axis#-- [m]
    #-- Ratio between ellipsoidal axes
    ratio = (1.0 -flat)
    #-- Polar radius of curvature
    pol_rad=a_axis/(1.0 -flat)

    #-- Linear eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    #-- first numerical eccentricity
    ecc1 = lin_ecc/a_axis
    #-- second numerical eccentricity
    ecc2 = lin_ecc/b_axis

    #-- m parameter [omega^2*a^2*b/(GM)]
    #-- p. 70, Eqn.(2-137)
    mp = omega**2*((1 -flat)*a_axis**3)/GM

    #-- q, q_0
    #-- p. 67, Eqn.(2-113)
    q = 0.5*((1.0 + 3.0/(ecc2**2))*np.arctan(ecc2)-3.0/ecc2)
    q_0 = 3*(1.0 +1.0/(ecc2**2))*(1.0 -1.0/ecc2*np.arctan(ecc2))-1.0

    #-- J_2 p. 75 Eqn.(2-167), p. 76 Eqn.(2-172)
    j_2 = (ecc1**2)*(1.0 - 2.0*mp*ecc2/(15.0*q))/3.0
    #-- Normalized C20 terms.
    #-- p. 60, Eqn.(2-80)
    C20 = -j_2/np.sqrt(5.0)

    #-- Normal gravity at the equator.
    #-- p. 71, Eqn.(2-141)
    ga = GM/(a_axis*b_axis)*(1.0 -mp -mp*ecc2*q_0/(6.0*q))

    #-- Normal gravity at the pole.
    #-- p. 71, Eqn.(2-142)
    gb = GM/(a_axis**2.0)*(1.0 +mp*ecc2*q_0/(3.0*q))

    #-- ratio between gravity at pole versus gravity at equator
    dk = b_axis*gb/(a_axis*ga) - 1.0

    #-- Normal potential at the ellipsoid
    #-- p. 68, Eqn.(2-123)
    U0 = GM/lin_ecc*np.arctan(ecc2)+(1.0/3.0)*omega**2*a_axis**2

    #-- Surface area of the reference ellipsoid [m^2]
    area = np.pi*a_axis**2.*(2.+((1.-ecc1**2)/ecc1)*np.log((1.+ecc1)/(1.-ecc1)))
    #-- Volume of the reference ellipsoid [m^3]
    vol = (4.0*np.pi/3.0)*(a_axis**3.0)*(1.0-ecc1**2.0)**0.5

    return {'a':a_axis, 'b':b_axis, 'f':flat, 'rad_p':pol_rad, 'rad_e':rad_e,
            'ratio':ratio, 'GM':GM, 'omega':omega, 'C20':C20, 'J2':j_2, 'U0':U0,
            'dk':dk, 'norm_a':ga, 'norm_b':gb, 'mp':mp, 'q':q, 'q0':q_0,
            'ecc':lin_ecc, 'ecc1':ecc1,'ecc2':ecc2, 'area':area, 'volume':vol}
