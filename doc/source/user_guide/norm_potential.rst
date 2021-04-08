=================
norm_potential.py
=================

- Calculates the normal potential at a given latitude and height
- Method described in [Barthelmes2013]_ and [Moazezi2012]_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.norm_potential import norm_potential
    U, dU_dr, dU_dtheta = norm_potential(lat, lon, h, 'WGS84', lmax)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/norm_potential.py

Arguments
#########

1. ``latitude``: latitude in degrees
2. ``longitude``: longitude in degrees
3. ``height``: height above reference ellipsoid in meters
4. ``refell``: reference ellipsoid name

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
5. ``lmax``: maximum spherical harmonic degree

Returns
#######

- ``U``: normal potential at height h
- ``dU_dr``: derivative of normal potential with respect to radius
- ``dU_dtheta``: derivative of normal potential with respect to theta

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [Moazezi2012] S. Moazezi and H. Zomorrodian, "GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly", *Earth Science Informatics*, 5, 123--136, (2012). `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_
