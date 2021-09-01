==========
spatial.py
==========

Utilities for operating on spatial data and converting ellipsoids

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/spatial.py

General Methods
===============

.. method:: geoid_toolkit.spatial.convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10)

    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Arguments:

        ``phi1``: latitude of input ellipsoid in degrees

        ``h1``: height above input ellipsoid in meters

        ``a1``: semi-major axis of input ellipsoid

        ``f1``: flattening of input ellipsoid

        ``a2``: semi-major axis of output ellipsoid

        ``f2``: flattening of output ellipsoid

    Keyword arguments:

        ``eps``: tolerance to prevent division by small numbers and to determine convergence

        ``itmax``: maximum number of iterations to use in Newton-Raphson

    Returns:

        ``phi2``: latitude of output ellipsoid in degrees

        ``h2``: height above output ellipsoid in meters


.. method:: geoid_toolkit.spatial.compute_delta_h(a1, f1, a2, f2, lat)

    Compute difference in elevation for two ellipsoids at a given latitude using a simplified empirical equation

    Arguments:

        ``a1``: semi-major axis of input ellipsoid

        ``f1``: flattening of input ellipsoid

        ``a2``: semi-major axis of output ellipsoid

        ``f2``: flattening of output ellipsoid

        ``lat``: array of latitudes in degrees

    Returns:

        ``delta_h``: difference in elevation for two ellipsoids

.. method:: geoid_toolkit.spatial.to_cartesian(lon,lat,a_axis=6378137.0,flat=1.0/298.257223563)

    Converts geodetic coordinates to Cartesian coordinates

    Arguments:

        ``lon``: longitude

        ``lat``: latitude

    Keyword arguments:

        ``h``: height

        ``a_axis``: semimajor axis of the ellipsoid

        ``flat``: ellipsoidal flattening

    Returns:

        ``x``, ``y``, ``z`` in Cartesian coordinates


.. method:: geoid_toolkit.spatial.to_sphere(x,y,z)

    Convert from Cartesian coordinates to spherical coordinates

    Arguments:

        ``x``, ``y``, ``z`` in Cartesian coordinates

    Returns:

        ``lon``: longitude

        ``lat``: latitude

        ``rad``: radius


.. method:: geoid_toolkit.spatial.to_geodetic(x,y,z,a_axis=6378137.0,flat=1.0/298.257223563)

    Convert from Cartesian coordinates to geodetic coordinates using `a closed form solution <https://arc.aiaa.org/doi/abs/10.2514/3.21016>`_

    Arguments:

        ``x``, ``y``, ``z`` in Cartesian coordinates

    Keyword arguments:

        ``a_axis``: semimajor axis of the ellipsoid

        ``flat``: ellipsoidal flattening

    Returns:

        ``lon``: longitude

        ``lat``: latitude

        ``h``: height
