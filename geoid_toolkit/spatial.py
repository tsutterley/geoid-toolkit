#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (09/2021)

Utilities for operating on spatial data and converting ellipsoids

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Written 09/2021
"""
import numpy as np

def convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10):
    """
    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Inputs:
        phi1: latitude of input ellipsoid in degrees
        h1: height above input ellipsoid in meters
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid

    Options:
        eps: tolerance to prevent division by small numbers
            and to determine convergence
        itmax: maximum number of iterations to use in Newton-Raphson

    Returns:
        phi2: latitude of output ellipsoid in degrees
        h2: height above output ellipsoid in meters

    References:
        Astronomical Algorithms, Jean Meeus, 1991, Willmann-Bell, Inc.
            pp. 77-82
    """
    if (len(phi1) != len(h1)):
        raise ValueError('phi and h have incompatable dimensions')
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- initialize output arrays
    npts = len(phi1)
    phi2 = np.zeros((npts))
    h2 = np.zeros((npts))
    #-- for each point
    for N in range(npts):
        #-- force phi1 into range -90 <= phi1 <= 90
        if (np.abs(phi1[N]) > 90.0):
            phi1[N] = np.sign(phi1[N])*90.0
        #-- handle special case near the equator
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + a1 - a2
        if (np.abs(phi1[N]) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + a1 - a2
        #-- handle special case near the poles
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + b1 - b2
        elif ((90.0 - np.abs(phi1[N])) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + b1 - b2
        #-- handle case if latitude is within 45 degrees of equator
        elif (np.abs(phi1[N]) <= 45):
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = b2 * b2 - a2 * a2
            k1 = a2 * hpr1cos
            k2 = b2 * hpr1sin
            #-- perform newton-raphson iteration to solve for u2
            #-- cos(u2) will not be close to zero since abs(phi1) <= 45
            for i in range(0, itmax+1):
                cosu2 = np.cos(u2)
                fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(phi2r)
        #-- handle final case where latitudes are between 45 degrees and pole
        else:
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = a2 * a2 - b2 * b2
            k1 = b2 * hpr1sin
            k2 = a2 * hpr1cos
            #-- perform newton-raphson iteration to solve for u2
            #-- sin(u2) will not be close to zero since abs(phi1) > 45
            for i in range(0, itmax+1):
                sinu2 = np.sin(u2)
                fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                fu2p =  -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(phi2r)

    #-- return the latitude and height
    return (phi2, h2)

def compute_delta_h(a1, f1, a2, f2, lat):
    """
    Compute difference in elevation for two ellipsoids at a given
        latitude using a simplified empirical equation

    Inputs:
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid
        lat: array of latitudes in degrees

    Returns:
        delta_h: difference in elevation for two ellipsoids

    Reference:
        J Meeus, Astronomical Algorithms, pp. 77-82 (1991)
    """
    #-- force phi into range -90 <= phi <= 90
    gt90, = np.nonzero((lat < -90.0) | (lat > 90.0))
    lat[gt90] = np.sign(lat[gt90])*90.0
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- compute delta_a and delta_b coefficients
    delta_a = a2 - a1
    delta_b = b2 - b1
    #-- compute differences between ellipsoids
    #-- delta_h = -(delta_a * cos(phi)^2 + delta_b * sin(phi)^2)
    phi = lat * np.pi/180.0
    delta_h = -(delta_a*np.cos(phi)**2 + delta_b*np.sin(phi)**2)
    return delta_h

def to_cartesian(lon,lat,h=0.0,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Converts geodetic coordinates to Cartesian coordinates

    Inputs:
        lon: longitude (degrees east)
        lat: latitude (degrees north)

    Options:
        h: height above ellipsoid (or sphere)
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
            * for spherical coordinates set to radius of the Earth
        flat: ellipsoidal flattening (default: WGS84)
            * for spherical coordinates set to 0
    """
    #-- verify axes
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    #-- fix coordinates to be 0:360
    count = np.count_nonzero(lon < 0)
    if (count != 0):
        lt0, = np.nonzero(lon < 0)
        lon[lt0] += 360.0
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- convert from geodetic latitude to geocentric latitude
    dtr = np.pi/180.0
    #-- geodetic latitude in radians
    latitude_geodetic_rad = lat*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(lon*dtr)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(lon*dtr)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    #-- return the cartesian coordinates
    return (X,Y,Z)

def to_sphere(x,y,z):
    """
    Convert from cartesian coordinates to spherical coordinates

    Inputs:
        x,y,z in cartesian coordinates
    """
    #-- calculate radius
    rad = np.sqrt(x**2.0 + y**2.0 + z**2.0)
    #-- calculate angular coordinates
    #-- phi: azimuthal angle
    phi = np.arctan2(y,x)
    #-- th: polar angle
    th = np.arccos(z/rad)
    #-- convert to degrees and fix to 0:360
    lon = 180.0*phi/np.pi
    count = np.count_nonzero(lon < 0)
    if (count != 0):
        lt0 = np.nonzero(lon < 0)
        lon[lt0] = lon[lt0]+360.0
    #-- convert to degrees and fix to -90:90
    lat = 90.0 - (180.0*th/np.pi)
    #-- return latitude, longitude and radius
    return (lon,lat,rad)

def to_geodetic(x,y,z,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using a closed form solution

    Inputs:
        x,y,z in cartesian coordinates

    Options:
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
        flat: ellipsoidal flattening (default: WGS84)

    References:
        J Zhu "Exact conversion of Earth-centered, Earth-fixed
            coordinates to geodetic coordinates"
        Journal of Guidance, Control, and Dynamics,
        16(2), 389--391, 1993
        https://arc.aiaa.org/doi/abs/10.2514/3.21016
    """
    #-- semiminor axis of the WGS84 ellipsoid [m]
    b_axis = (1.0 - flat)*a_axis
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- square of first numerical eccentricity
    e12 = ecc1**2
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- calculate distance
    w = np.sqrt(x**2 + y**2)
    #-- calculate longitude
    lon = np.arctan2(y,x)/dtr
    lat = np.zeros_like(lon)
    h = np.zeros_like(lon)
    if (w == 0):
        #-- special case where w == 0 (exact polar solution)
        h = np.sign(z)*z - b_axis
        lat = 90.0*np.sign(z)
    else:
        #-- all other cases
        l = e12/2.0
        m = (w/a_axis)**2.0
        n = ((1.0-e12)*z/b_axis)**2.0
        i = -(2.0*l**2 + m + n)/2.0
        k = (l**2.0 - m - n)*l**2.0
        q = (1.0/216.0)*(m + n - 4.0*l**2)**3.0 + m*n*l**2.0
        D = np.sqrt((2.0*q - m*n*l**2)*m*n*l**2)
        B = i/3.0 - (q+D)**(1.0/3.0) - (q-D)**(1.0/3.0)
        t = np.sqrt(np.sqrt(B**2-k) - (B+i)/2.0)-np.sign(m-n)*np.sqrt((B-i)/2.0)
        wi = w/(t+l)
        zi = (1.0-e12)*z/(t-l)
        #-- calculate latitude and height
        lat = np.arctan2(zi,((1.0-e12)*wi))/dtr
        h = np.sign(t-1.0+l)*np.sqrt((w-wi)**2.0 + (z-zi)**2.0)
    #-- return latitude, longitude and height
    return (lon,lat,h)
