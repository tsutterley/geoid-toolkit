=========================
calculate_tidal_offset.py
=========================

- Calculates the spherical harmonic offset for a tide system to change from a tide free state where there is no permanent direct and indirect tidal potentials [Losch2003]_
- Method described in [Barthelmes2013]_ and [Moazezi2012]_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.calculate_tidal_offset import calculate_tidal_offset
    delta = calculate_tidal_offset(TIDE, GM, R, refell)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/calculate_tidal_offset.py

Arguments
#########

1. ``TIDE``: output tidal system

    * ``'mean_tide'``: restores permanent tidal potentials (direct and indirect)
    * ``'zero_tide'``: restores permanent direct tidal potential
2. ``R``: average radius used in gravity model
3. ``GM``: geocentric graviational constant used in gravity model
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

Returns
#######

- ``delta``: offset for changing from tide free system

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [Losch2003] M. Losch and V. Seufer, "How to Compute Geoid Undulations (Geoid Height Relative to a Given Reference Ellipsoid) from Spherical Harmonic Coefficients for Satellite Altimetry Applications", (2003). `eprint ID: 11802 <http://mitgcm.org/~mlosch/geoidcookbook.pdf>`_

.. [Moazezi2012] S. Moazezi and H. Zomorrodian, "GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly", *Earth Science Informatics*, 5, 123--136, (2012). `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_
