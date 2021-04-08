================
ref_ellipsoid.py
================

- Computes parameters for a reference ellipsoid [HofmannWellenhof2006]_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.ref_ellipsoid import ref_ellipsoid
    params = ref_ellipsoid(refell)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/ref_ellipsoid.py

Arguments
#########

 1. ``refell``: reference ellipsoid name

    * ``'CLK66'``: Clarke 1866
    * ``'GRS67'``: Geodetic Reference System 1967
    * ``'GRS80'``: Geodetic Reference System 1980
    * ``'WGS72'``: World Geodetic System 1972
    * ``'WGS84'``: World Geodetic System 1984
    * ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
    * ``'NAD27'``: North American Datum 1927 (=CLK66)
    * ``'NAD83'``: North American Datum 1983 (=GRS80)
    * ``'INTER'``: International
    * ``'KRASS'``: Krassovsky (USSR)
    * ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
    * ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
    * ``'EGM96'``: EGM 1996 gravity model

Keyword arguments
#################

- ``UNITS``: output units
    * ``'MKS'``: meters, kilograms, seconds
    * ``'CGS'``: centimeters, grams, seconds

Returns
#######

- ``a``: semimajor semi-axis (m)
- ``b``: semiminor semi-axis (m)
- ``f`` : flattening
- ``c`` : Polar radius of curvature
- ``rad_e``: mean radius of ellipsoid having the same volume
- ``rad_p``: Polar radius of curvature
- ``C20``: Normalized C20 harmonic
- ``norm_a``: Normal gravity at the equator
- ``norm_b``: Normal gravity at the pole
- ``U0``: Normal potential at the ellipsoid
- ``dk``: ratio between gravity at pole versus gravity at equator
- ``m``: m parameter (m)
- ``lin_ecc``: Linear eccentricity
- ``ecc1``: First eccentricity
- ``ecc2``: Second eccentricity
- ``area``: Area of the ellipsoid

References
##########

.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
