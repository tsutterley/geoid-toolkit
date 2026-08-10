#!/usr/bin/env python
"""
ref_ellipsoid.py
Written by Tyler Sutterley (12/2022)

Computes parameters for a reference ellipsoid

CALLING SEQUENCE
    wgs84 = ref_ellipsoid('WGS84', UNITS='MKS')

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
    rad_e: mean radius of ellipsoid having the same volume
    rad_p: Polar radius of curvature
    C20: Normalized C20 harmonic
    J2: Oblateness coefficient
    norm_a: Normal gravity at the equator
    norm_b: Normal gravity at the pole
    U0: Normal potential at the ellipsoid
    dk: ratio between gravity at pole versus gravity at equator
    m: m parameter (m)
    lin_ecc: Linear eccentricity
    ecc1: First eccentricity
    ecc2: Second eccentricity
    area: Area of the ellipsoid
    volume: Volume of the ellipsoid
    rho_e: Average density

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Hofmann-Wellenhof and Moritz (2006)

UPDATE HISTORY:
    Updated 12/2022: output average Earth's density
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

from __future__ import annotations

import numpy as np
import geoid_toolkit.math


__all__ = [
    'ref_ellipsoid',
    'norm_gravity',
    'norm_potential',
]


def ref_ellipsoid(refell: str, UNITS: str = 'MKS') -> dict:
    """
    Computes parameters for a reference ellipsoid
    :cite:p:`HofmannWellenhof:2006hy`

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
    rad_e: float
        Mean radius of ellipsoid having the same volume
    rad_p: float
        Polar radius of curvature
    C20: float
        Normalized C20 harmonic
    J2: float
        Oblateness coefficient
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
    volume: float
        Volume of the ellipsoid
    rho_e: float
        Average density
    """
    # validate units
    assert UNITS in ('MKS', 'CGS')

    # set parameters for ellipsoid
    if refell.upper() in ('CLK66', 'NAD27'):
        # Clarke 1866
        a_axis = 6378206.4  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 294.9786982  # flattening of the ellipsoid

    elif refell.upper() in ('GRS80', 'NAD83'):
        # Geodetic Reference System 1980
        # North American Datum 1983
        a_axis = 6378135.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.26  # flattening of the ellipsoid
        GM = 3.986005e14  # [m^3/s^2] Geocentric Gravitational Constant

    elif refell.upper() == 'GRS67':
        # Geodetic Reference System 1967
        # International Astronomical Union (IAU ellipsoid)
        a_axis = 6378160.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.247167427  # flattening of the ellipsoid
        GM = 3.98603e14  # [m^3/s^2] Geocentric Gravitational Constant
        omega = 7292115.1467e-11  # angular velocity of the Earth [rad/s]

    elif refell.upper() == 'WGS72':
        # World Geodetic System 1972
        a_axis = 6378135.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.26  # flattening of the ellipsoid

    elif refell.upper() == 'WGS84':
        # World Geodetic System 1984
        a_axis = 6378137.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.257223563  # flattening of the ellipsoid

    elif refell.upper() == 'ATS77':
        # Quasi-earth centred ellipsoid for ATS77
        a_axis = 6378135.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.257  # flattening of the ellipsoid

    elif refell.upper() == 'KRASS':
        # Krassovsky (USSR)
        a_axis = 6378245.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.3  # flattening of the ellipsoid

    elif refell.upper() == 'INTER':
        # International
        a_axis = 6378388.0  # [m] semimajor axis of the ellipsoid
        flat = 1 / 297.0  # flattening of the ellipsoid

    elif refell.upper() == 'MAIRY':
        # Modified Airy (Ireland 1965/1975)
        a_axis = 6377340.189  # [m] semimajor axis of the ellipsoid
        flat = 1 / 299.3249646  # flattening of the ellipsoid

    elif refell.upper() == 'TOPEX':
        # TOPEX/POSEIDON ellipsoid
        a_axis = 6378136.3  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.257  # flattening of the ellipsoid
        GM = 3.986004415e14  # [m^3/s^2]

    elif refell.upper() == 'EGM96':
        # EGM 1996 gravity model
        a_axis = 6378136.3  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.256415099  # flattening of the ellipsoid
        GM = 3.986004415e14  # [m^3/s^2]

    elif refell.upper() == 'HGH80':
        # Hughes 1980 Ellipsoid used in some NSIDC data
        a_axis = 6378273.0  # [m] semimajor axis of the ellipsoid
        flat = 1.0 / 298.279411123064  # flattening of the ellipsoid

    else:
        raise ValueError('Incorrect reference ellipsoid Name')

    if refell.upper() not in ('GRS80', 'GRS67', 'NAD83', 'TOPEX', 'EGM96'):
        # for ellipsoids not listing the Geocentric Gravitational Constant
        GM = 3.986004418e14  # [m^3/s^2]

    if refell.upper() not in ('GRS67'):
        # for ellipsoids not listing the angular velocity of the Earth
        omega = 7292115e-11  # [rad/s]

    # universal gravitational constant [N*m^2/kg^2]
    G = 6.67430e-11

    # convert units to CGS
    if UNITS == 'CGS':
        a_axis *= 100.0
        GM *= 1e6
        G *= 1000.0  # [dyn*cm^2/g^2]

    # DERIVED PARAMETERS:
    # mean radius of the Earth having the same volume
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 - f)
    rad_e = a_axis * (1.0 - flat) ** (1.0 / 3.0)

    # semiminor axis of the ellipsoid
    b_axis = (1.0 - flat) * a_axis  # [m]
    # Ratio between ellipsoidal axes
    ratio = 1.0 - flat
    # Polar radius of curvature
    pol_rad = a_axis / (1.0 - flat)

    # Linear eccentricity
    lin_ecc = np.sqrt((2.0 * flat - flat**2) * a_axis**2)
    # first numerical eccentricity
    ecc1 = lin_ecc / a_axis
    # second numerical eccentricity
    ecc2 = lin_ecc / b_axis

    # m parameter [omega^2*a^2*b/(GM)]
    # p. 70, Eqn.(2-137)
    mp = omega**2 * ((1 - flat) * a_axis**3) / GM

    # q, q_0
    # p. 67, Eqn.(2-113)
    q = 0.5 * ((1.0 + 3.0 / (ecc2**2)) * np.arctan(ecc2) - 3.0 / ecc2)
    q_0 = (
        3 * (1.0 + 1.0 / (ecc2**2)) * (1.0 - 1.0 / ecc2 * np.arctan(ecc2)) - 1.0
    )

    # J_2 p. 75 Eqn.(2-167), p. 76 Eqn.(2-172)
    j_2 = (ecc1**2) * (1.0 - 2.0 * mp * ecc2 / (15.0 * q)) / 3.0
    # Normalized C20 terms.
    # p. 60, Eqn.(2-80)
    C20 = -j_2 / np.sqrt(5.0)

    # Normal gravity at the equator.
    # p. 71, Eqn.(2-141)
    ga = GM / (a_axis * b_axis) * (1.0 - mp - mp * ecc2 * q_0 / (6.0 * q))

    # Normal gravity at the pole.
    # p. 71, Eqn.(2-142)
    gb = GM / (a_axis**2.0) * (1.0 + mp * ecc2 * q_0 / (3.0 * q))

    # ratio between gravity at pole versus gravity at equator
    dk = b_axis * gb / (a_axis * ga) - 1.0

    # Normal potential at the ellipsoid
    # p. 68, Eqn.(2-123)
    U0 = GM / lin_ecc * np.arctan(ecc2) + (1.0 / 3.0) * omega**2 * a_axis**2

    # Surface area of the reference ellipsoid [m^2]
    area = (
        np.pi
        * a_axis**2.0
        * (2.0 + ((1.0 - ecc1**2) / ecc1) * np.log((1.0 + ecc1) / (1.0 - ecc1)))
    )
    # Volume of the reference ellipsoid [m^3]
    vol = (4.0 * np.pi / 3.0) * (a_axis**3.0) * (1.0 - ecc1**2.0) ** 0.5
    # Average density [kg/m^3]
    rho_e = GM / (G * vol)

    return {
        'a': a_axis,
        'b': b_axis,
        'f': flat,
        'rad_p': pol_rad,
        'rad_e': rad_e,
        'ratio': ratio,
        'GM': GM,
        'omega': omega,
        'C20': C20,
        'J2': j_2,
        'U0': U0,
        'dk': dk,
        'norm_a': ga,
        'norm_b': gb,
        'mp': mp,
        'q': q,
        'q0': q_0,
        'ecc': lin_ecc,
        'ecc1': ecc1,
        'ecc2': ecc2,
        'area': area,
        'volume': vol,
        'rho_e': rho_e,
    }


def norm_gravity(lat: float | np.ndarray, h: float | np.ndarray, refell: str):
    """
    Calculates the normal gravity of an ellipsoid and calculates the derivative
    with respect to height following :cite:t:`HofmannWellenhof:2006hy`

    Parameters
    ----------
    lat: float
        latitude in degrees
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

    Returns
    -------
    gamma_h: float
        normal gravity for ellipsoid at height
    dgamma_dh: float
        derivative of normal gravity with respect to height
    """

    # convert latitude from degrees to radians
    phi = np.radians(lat)

    # get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = ellip['a']
    b = ellip['b']
    # eccentricity
    ecc2 = ellip['ecc2']
    GM = ellip['GM']
    # m parameter [omega^2*a^2*b/(GM)]
    m = ellip['mp']
    # flattening components
    f = ellip['f']
    f_2 = (
        -f
        + (5.0 / 2.0) * m
        + (1.0 / 2.0) * f**2.0
        - (26.0 / 7.0) * f * m
        + (15.0 / 4.0) * m**2.0
    )
    f_4 = -(1.0 / 2.0) * f**2.0 + (5.0 / 2.0) * f * m

    # Normal gravity at the equator.
    # p. 79, Eqn.(2-186)
    gamma_a = (GM / (a * b)) * (
        1.0 - (3.0 / 2.0) * m - (3.0 / 14.0) * ecc2**2.0 * m
    )
    # Normal gravity
    # p. 80, Eqn.(2-199)
    gamma_0 = gamma_a * (
        1.0 + f_2 * np.sin(phi) ** 2.0 + f_4 * np.sin(phi) ** 4.0
    )
    # Normal gravity at height
    # p. 82, Eqn.(2-215)
    p_1 = 1.0 + f + m - 2.0 * f * np.sin(phi) ** 2.0
    gamma_h = gamma_0 * (1.0 - (2.0 / a) * p_1 * h + (3.0 / (a**2.0)) * h**2.0)
    # approximate derivative of normal gravity with respect to height
    dgamma_dh = ((-2.0 * gamma_0) / a) * p_1
    # return the normal gravity and the derivative
    return (gamma_h, dgamma_dh)


def norm_potential(
    lat: float | np.ndarray,
    lon: float | np.ndarray,
    h: float | np.ndarray,
    refell: str,
    lmax: int,
):
    """
    Calculates the normal potential following
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
    lmax: int
        maximum spherical harmonic degree

    Returns
    -------
    U: float
        normal potential at height
    dU_dr: float
        derivative of normal potential with respect to radius
    dU_dtheta: float
        derivative of normal potential with respect to theta
    """
    # import function to convert from geodetic to cartesian coordinates
    from geoid_toolkit.spatial import to_cartesian

    # get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = np.longdouble(ellip['a'])
    ecc1 = np.longdouble(ellip['ecc1'])
    GM = np.longdouble(ellip['GM'])
    J2 = np.longdouble(ellip['J2'])

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

    # calculate even zonal harmonics
    n = np.arange(2, 12 + 2, 2, dtype=np.longdouble)
    J2n = _cosine_even_zonals(J2, ecc1, n / 2.0)
    # normalized cosine harmonics: Cn = -Jn/np.sqrt(2.0*n+1.0)
    # J2 = 0.108262982131e-2
    C_2 = -J2n[0] / np.sqrt(5.0)
    # J4 = -0.237091120053e-5
    C_4 = -J2n[1] / np.sqrt(9.0)
    # J6 = 0.608346498882e-8
    C_6 = -J2n[2] / np.sqrt(13.0)
    # J8 = -0.142681087920e-10
    C_8 = -J2n[3] / np.sqrt(17.0)
    # J10 = 0.121439275882e-13
    C_10 = -J2n[4] / np.sqrt(21.0)
    # J12 = 0.205395070709e-15
    C_12 = -J2n[5] / np.sqrt(25.0)

    # calculate legendre polynomials at latitude and their first derivative
    Pl, dPl = geoid_toolkit.math.legendre_polynomials(
        lmax, np.cos(theta), ASTYPE=np.longdouble
    )

    # normal potentials and derivatives
    U = (GM / rr) * (
        1.0
        + (a / rr) ** 2.0 * C_2 * Pl[2, :]
        + (a / rr) ** 4.0 * C_4 * Pl[4, :]
        + (a / rr) ** 6.0 * C_6 * Pl[6, :]
        + (a / rr) ** 8.0 * C_8 * Pl[8, :]
        + (a / rr) ** 10.0 * C_10 * Pl[10, :]
        + (a / rr) ** 12.0 * C_12 * Pl[12, :]
    )
    dU_dr = GM * (
        -1.0 / rr**2.0
        - 3.0 * (a**2.0 / rr**4.0) * C_2 * Pl[2, :]
        - 5.0 * (a**4.0 / rr**6.0) * C_4 * Pl[4, :]
        - 7.0 * (a**6.0 / rr**8.0) * C_6 * Pl[6, :]
        - 9.0 * (a**8.0 / rr**10.0) * C_8 * Pl[8, :]
        - 11.0 * (a**10.0 / rr**12.0) * C_10 * Pl[10, :]
        - 13.0 * (a**12.0 / rr**14.0) * C_12 * Pl[12, :]
    )
    dU_dtheta = (GM / rr) * (
        1.0
        + (a / rr) ** 2.0 * C_2 * dPl[2, :]
        + (a / rr) ** 4.0 * C_4 * dPl[4, :]
        + (a / rr) ** 6.0 * C_6 * dPl[6, :]
        + (a / rr) ** 8.0 * C_8 * dPl[8, :]
        + (a / rr) ** 10.0 * C_10 * dPl[10, :]
        + (a / rr) ** 12.0 * C_12 * dPl[12, :]
    )

    # return the potentials
    return (U, dU_dr, dU_dtheta)


# PURPOSE: Calculate even zonal harmonics using J2 and first eccentricity
def _cosine_even_zonals(
    J2: float,
    e: float,
    n: int,
) -> float:
    """
    Calculate even zonal harmonics using J2 and first eccentricity

    Parameters
    ----------
    J2: float
        Oblateness
    e: float
        First eccentricity
    n: int
        spherical harmonic degree

    Returns
    -------
    J2n: float
        Even zonal harmonics
    """
    # p. 76 Eqn.(2-170)
    J2n = (
        (-1.0) ** (n + 1.0)
        * ((3.0 * e ** (2.0 * n)) / ((2.0 * n + 1.0) * (2.0 * n + 3.0)))
        * (1.0 - n + 5.0 * n * J2 / (e**2.0))
    )
    return J2n
