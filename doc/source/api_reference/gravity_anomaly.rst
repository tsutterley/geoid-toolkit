===============
gravity_anomaly
===============

- Calculates the gravity anomaly at a given latitude and longitude using different methods

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.gravity_anomaly import gravity_anomaly
    ddelta_g = gravity_anomaly(lat, lon, h, 'WGS84', clm, slm, lmax, R, GM, METHOD='first')

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/gravity_anomaly.py

.. autofunction:: geoid_toolkit.gravity_anomaly
