================
geoid_undulation
================

- Calculates the geoidal undulation at a given latitude and longitude using an iterative approach

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.geoid_undulation import geoid_undulation
    N = geoid_undulation(lat, lon, 'WGS84', clm, slm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/geoid_undulation.py

.. autofunction:: geoid_toolkit.geoid_undulation
