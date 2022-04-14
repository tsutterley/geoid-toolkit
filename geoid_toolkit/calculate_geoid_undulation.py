#!/usr/bin/env python
u"""
calculate_geoid_undulation.py
Written by Tyler Sutterley (04/2022)
Wrapper function for computing geoid undulations from a gravity model

INPUTS:
    lon: longitudinal points to calculate geoid height
    lat: latitudinal points to calculate geoid height
    gravity_model_file: full path to static gravity model file

OPTIONS:
    ELLIPSOID: reference ellipsoid name
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
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data
    LMAX: maximum spherical harmonic degree (level of truncation)
    TIDE: tide system of output geoid
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
            this is the default (leaving the model as is)
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)
    EPS: level of precision for calculating geoid height
    ZIP: input gravity field file is compressed in an archive file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    geoid_undulation.py: geoidal undulation at a given latitude and longitude
    read_ICGEM_harmonics.py: reads the coefficients for a given gravity model file
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: add more keyword options to match read ICGEM options
    Updated 09/2021: define int/float precision to prevent deprecation warning
    Updated 11/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.geoid_undulation import geoid_undulation
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics

#-- PURPOSE: calculate geoid heights at a set of latitudes and longitudes
def calculate_geoid_undulation(lon, lat, gravity_model_file, **kwargs):
    """
    Wrapper function for computing geoid undulations from a gravity model

    Parameters
    ----------
    lon: float
        longitudinal points to calculate geoid height
    lat: float
        latitudinal points to calculate geoid height
    gravity_model_file: str
        full path to static gravity model file
    LMAX: int or NoneType, default None
        maximum spherical harmonic degree
    ELLIPSOID: str, default 'WGS84'
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
    TIDE: str, default 'tide_free'
        Permanent tide system of output geoid

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    GAUSS: int, default 0
        Gaussian Smoothing Radius in km
    EPS: float, default 1e-8
        level of precision for calculating geoid height
    ZIP: bool, default False
        Gravity field file is compressed in an archive file

    Returns
    -------
    N: float
        geoidal undulation for a given ellipsoid in meters
    """
    #-- set default keyword arguments
    kwargs.setdefault('LMAX',None)
    kwargs.setdefault('ELLIPSOID','WGS84')
    kwargs.setdefault('TIDE','tide_free')
    kwargs.setdefault('GAUSS',0)
    kwargs.setdefault('EPS',1e-8)
    kwargs.setdefault('ZIP',False)
    #-- read gravity model Ylms and change tide if specified
    Ylms = read_ICGEM_harmonics(gravity_model_file,**kwargs)
    R = np.float64(Ylms['radius'])
    GM = np.float64(Ylms['earth_gravity_constant'])
    LMAX = np.int64(Ylms['max_degree'])
    #-- calculate geoid at coordinates
    N = geoid_undulation(lat, lon, kwargs['ELLIPSOID'],
        Ylms['clm'], Ylms['slm'], LMAX, R, GM,
        GAUSS=kwargs['GAUSS'], EPS=kwargs['EPS'])
    #-- return the geoid undulation
    return N
