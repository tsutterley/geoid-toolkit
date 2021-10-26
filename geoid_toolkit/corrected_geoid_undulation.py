#!/usr/bin/env python
u"""
corrected_geoid_undulation.py
Written by Tyler Sutterley (07/2017)
Calculates the topographically corrected geoidal undulation at a given latitude
    and longitude using an iterative approach described in Barthelmes (2009)

CALLING SEQUENCE:
    geoid = corrected_geoid_undulation(lat, lon, 'WGS84', clm, slm, tclm, tslm,
        lmax, R, GM, density)

INPUT:
    latitude: latitude in degrees
    longitude: latitude in degrees
    refell: reference ellipsoid specified in ref_ellipsoid.py
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
    clm: cosine spherical harmonics for a gravity model
    slm: sine spherical harmonics for a gravity model
    tclm: cosine spherical harmonics for a topographic model
    tslm: sine spherical harmonics for a topographic model
    lmax: maximum spherical harmonic degree
    R: average radius used in gravity model
    GM: geocentric graviational constant used in gravity model
    density: density of the topography in the model

OUTPUT:
    geoidal undulation for a given ellipsoid in meters

OPTIONS:
    GAUSS: Gaussian Smoothing Radius in km (default is no filtering)
    EPS: level of precision for calculating geoid height

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    real_potential.py: real potential at lat, lon and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    topographic_potential.py: topographic potential correction at lat, lon
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

REFERENCES:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4
    Barthelmes, "Definition of Functionals of the Geopotential and Their
        Calculation from Spherical Harmonic Models", STR09/02 (2009)
        http://icgem.gfz-potsdam.de/ICGEM/theory/str-0902-revised.pdf

UPDATE HISTORY:
    Written 07/2017
"""
import numpy as np
from geoid_toolkit.real_potential import real_potential
from geoid_toolkit.norm_potential import norm_potential
from geoid_toolkit.topographic_potential import topographic_potential
from geoid_toolkit.norm_gravity import norm_gravity

def corrected_geoid_undulation(lat, lon, refell, clm, slm, tclm, tslm, lmax,
    R, GM, density, GAUSS=0, EPS=1e-8):
    #-- calculate the real and normal potentials for the first iteration
    W,dWdr=real_potential(lat,lon,0.0,refell,GM,R,clm,slm,lmax,GAUSS=GAUSS)
    U,dUdr,dUdt = norm_potential(lat,lon,0.0,refell,lmax)
    #-- topographic potential correction
    T = topographic_potential(lon,lat,refell,R,tclm,tslm,lmax,density,GAUSS=GAUSS)
    #-- normal gravity at latitude
    gamma_h,dgamma_dh = norm_gravity(lat, 0.0, refell)
    #-- geoid height for first iteration
    N_1 = (W - U - T) / gamma_h
    #-- set geoid height to the first iteration and set RMS as infinite
    N = np.copy(N_1)
    RMS = np.inf
    while (RMS > EPS):
        #-- calculate the real potentials for the iteration
        W,dWdr = real_potential(lat,lon,N_1,refell,GM,R,clm,slm,lmax, GAUSS=GAUSS)
        #-- add geoid height for iteration
        N_1 += (W - U - T) / gamma_h
        #-- calculate RMS between iterations
        RMS = np.sqrt(np.sum((N - N_1)**2)/len(lat))
        #-- set N to the previous iteration
        N = np.copy(N_1)
    #-- return the geoid height
    return N
