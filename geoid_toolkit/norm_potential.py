#!/usr/bin/env python
u"""
norm_potential.py
Written by Tyler Sutterley (04/2022)
Calculates the normal potential at a given latitude and height

CALLING SEQUENCE:
    U, dU_dr, dU_dtheta = norm_potential(lat, lon, h, 'WGS84', lmax)

INPUT:
    lat: latitude in degrees
    lon: longitude in degrees
    h: height above reference ellipsoid in meters
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
    lmax: maximum spherical harmonic degree

OUTPUT:
    U: normal potential at height h
    dU_dr: derivative of normal potential with respect to radius
    dU_dtheta: derivative of normal potential with respect to theta

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    legendre_polynomials.py: Computes fully normalized Legendre polynomials

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
    Updated 07/2017: changed dtypes to long double for high degree models
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
from geoid_toolkit.legendre_polynomials import legendre_polynomials

def norm_potential(lat, lon, h, refell, lmax):
    """
    Calculates the normal potential at a given latitude and height
    following [Barthelmes2013]_ and [Moazezi2012]_

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

    #-- get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    a = np.longdouble(ellip['a'])
    ecc1 = np.longdouble(ellip['ecc1'])
    GM = np.longdouble(ellip['GM'])
    J2 = np.longdouble(ellip['J2'])

    #-- convert from geodetic latitude to geocentric latitude
    latitude_geodetic_rad = (np.pi*lat/180.0).astype(np.longdouble)
    longitude_rad = (np.pi*lon/180.0).astype(np.longdouble)
    N = a/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(longitude_rad)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))

    #-- calculate even zonal harmonics
    n = np.arange(2, 12+2, 2, dtype=np.longdouble)
    J2n = cosine_even_zonals(J2, ecc1, n/2.0)
    #-- normalized cosine harmonics: Cn = -Jn/np.sqrt(2.0*n+1.0)
    #-- J2 = 0.108262982131e-2
    C_2 = -J2n[0]/np.sqrt(5.0)
    #-- J4 = -0.237091120053e-5
    C_4 = -J2n[1]/np.sqrt(9.0)
    #-- J6 = 0.608346498882e-8
    C_6 = -J2n[2]/np.sqrt(13.0)
    #-- J8 = -0.142681087920e-10
    C_8 = -J2n[3]/np.sqrt(17.0)
    #-- J10 = 0.121439275882e-13
    C_10 = -J2n[4]/np.sqrt(21.0)
    #-- J12 = 0.205395070709e-15
    C_12 = -J2n[5]/np.sqrt(25.0)

    #-- calculate legendre polynomials at latitude and their first derivative
    Pl,dPl = legendre_polynomials(lmax, np.sin(latitude_geocentric),
        ASTYPE=np.longdouble)

    #-- normal potentials and derivatives
    U = (GM/rr) * (1.0 + (a/rr)**2.*C_2*Pl[2,:] + (a/rr)**4.*C_4*Pl[4,:] + \
        (a/rr)**6.*C_6*Pl[6,:] + (a/rr)**8.*C_8*Pl[8,:] + \
        (a/rr)**10.*C_10*Pl[10,:] + (a/rr)**12.*C_12*Pl[12,:])
    dU_dr = GM * (-1.0 / rr**2.0 - 3.0*(a**2.0/rr**4.0)*C_2*Pl[2,:] - \
        5.0*(a**4.0/rr**6.0)*C_4*Pl[4,:] -7.0*(a**6.0/rr**8.0)*C_6*Pl[6,:] - \
        9.0*(a**8.0/rr**10.)*C_8*Pl[8,:] -11.*(a**10./rr**12.)*C_10*Pl[10,:] - \
        13.*(a**12./rr**14.)*C_12*Pl[12,:])
    dU_dtheta = (GM/rr) * (1.0 + (a/rr)**2.0*C_2*dPl[2,:] + \
        (a/rr)**4.0*C_4*dPl[4,:] + (a/rr)**6.0*C_6*dPl[6,:] + \
        (a/rr)**8.0*C_8*dPl[8,:] + (a/rr)**10.0*C_10*dPl[10,:] + \
        (a/rr)**12.0*C_12*dPl[12,:])

    #-- return the potentials
    return (U, dU_dr, dU_dtheta)

#-- PURPOSE: Calculate even zonal harmonics using J2 and first eccentricity
def cosine_even_zonals(J2,e,n):
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
    #-- p. 76 Eqn.(2-170)
    J2n = (-1.0)**(n+1.0)*((3.0*e**(2.0*n))/((2.0*n + 1.0)*(2.0*n + 3.0))) * \
        (1.0 - n + 5.0*n*J2/(e**2.0))
    return J2n
