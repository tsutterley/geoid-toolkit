=================
height_anomaly.py
=================

- Calculates the height anomaly at a given latitude and longitude using an iterative approach

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.height_anomaly import height_anomaly
    zeta = height_anomaly(lat, lon, h, 'WGS84', clm, slm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/height_anomaly.py

.. autofunction:: geoid_toolkit.height_anomaly
