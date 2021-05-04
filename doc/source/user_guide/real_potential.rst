=================
real_potential.py
=================

- Calculates the real potential at a given latitude and height using coefficients from a gravity model
- Method described in [Barthelmes2013]_ and [Moazezi2012]_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.real_potential import real_potential
    W, dW_dr, dW_dtheta = real_potential(lat, lon, h, clm, slm, lmax)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/real_potential.py

Arguments
#########

1. ``latitude``: latitude in degrees
2. ``longitude``: longitude in degrees
3. ``height``: height above reference ellipsoid in meters
4. ``GM``: geocentric graviational constant used in gravity model
5. ``R``: average radius used in gravity model
6. ``clm``: cosine spherical harmonics for a gravity model
7. ``slm``: sine spherical harmonics for a gravity model
8. ``lmax``: maximum spherical harmonic degree

Keyword arguments
#################

- ``GAUSS``: Gaussian Smoothing Radius in km

Returns
#######

- ``W``: real potential at height h
- ``dW_dr``: derivative of real potential with respect to radius
- ``dW_dtheta``: derivative of real potential with respect to theta

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [Moazezi2012] S. Moazezi and H. Zomorrodian, "GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly", *Earth Science Informatics*, 5, 123--136, (2012). `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_
