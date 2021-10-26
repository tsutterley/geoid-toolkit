=============================
corrected_geoid_undulation.py
=============================

- Calculates the topographically corrected geoidal undulation at a given latitude and longitude using an iterative approach
- Method described in [Barthelmes2013]_ and [Moazezi2012]_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.corrected_geoid_undulation import corrected_geoid_undulation
    N = corrected_geoid_undulation(lat, lon, 'WGS84', clm, slm, tclm, tslm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/corrected_geoid_undulation.py

Arguments
#########

1. ``lat``: latitudinal points to calculate geoid height
2. ``lon``: longitudinal points to calculate geoid height
3. ``refell``: reference ellipsoid name

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
4. ``clm``: cosine spherical harmonics for a gravity model
5. ``slm``: sine spherical harmonics for a gravity model
6. ``tclm``: cosine spherical harmonics for a topographic model
7. ``tslm``: sine spherical harmonics for a topographic model
8. ``lmax``: maximum spherical harmonic degree
9. ``R``: average radius used in gravity model
10. ``GM``: geocentric graviational constant used in gravity model
11. ``density``: density of the topography in the model

Keyword arguments
#################

- ``GAUSS``: Gaussian Smoothing Radius in km
- ``EPS``: level of precision for calculating geoid height

Returns
#######

- ``N``: geoidal undulation for a given ellipsoid in meters

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [Moazezi2012] S. Moazezi and H. Zomorrodian, "GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly", *Earth Science Informatics*, 5, 123--136, (2012). `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_
