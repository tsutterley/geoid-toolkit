#!/usr/bin/env python
u"""
calculate_geoid_undulation.py
Written by Tyler Sutterley (11/2020)
Wrapper function for computing geoid undulations from a gravity model

INPUTS:
    lon: longitudinal points to calculate geoid height
    lat: latitudinal points to calculate geoid height
    gravity_model_file: full path to static gravity model file

OPTIONS:
    REFERENCE: reference ellipsoid name
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
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data
    LMAX: maximum spherical harmonic degree (level of truncation)
    TIDE: tide system of output geoid
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
            this is the default (leaving the model as is)
        mean_tide: restores permanent tidal potentials (direct and indirect)
        zero_tide: restores permanent direct tidal potential
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    geoid_undulation.py: geoidal undulation at a given latitude and longitude
    read_gravity_model.py: reads the coefficients for a given gravity model file
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

UPDATE HISTORY:
    Updated 11/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.geoid_undulation import geoid_undulation
from geoid_toolkit.read_gravity_model import read_gravity_model

#-- PURPOSE: calculate geoid heights at a set of latitudes and longitudes
def calculate_geoid_undulation(lon, lat, gravity_model_file, REFERENCE='WGS84',
    LMAX=None, TIDE='tide_free', GAUSS=0):
    """
    Wrapper function for computing geoid undulations from a gravity model

    Arguments
    ---------
    lon: longitudinal points to calculate geoid height
    lat: latitudinal points to calculate geoid height
    gravity_model_file: full path to static gravity model file

    Keyword arguments
    -----------------
    REFERENCE: reference ellipsoid name
    LMAX: maximum spherical harmonic degree (level of truncation)
    TIDE: tide system of output geoid
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)

    Returns
    -------
    N: geoidal undulation for a given ellipsoid in meters
    """
    #-- read gravity model Ylms and change tide if specified
    Ylms = read_gravity_model(gravity_model_file,LMAX,TIDE=TIDE)
    R = np.float(Ylms['radius'])
    GM = np.float(Ylms['earth_gravity_constant'])
    LMAX = np.int(Ylms['max_degree']) if not LMAX else LMAX
    #-- calculate geoid at coordinates
    N = geoid_undulation(lat, lon, REFERENCE, Ylms['clm'], Ylms['slm'], LMAX,
        R, GM, GAUSS=GAUSS, EPS=1e-8)
    #-- return the geoid undulation
    return N
