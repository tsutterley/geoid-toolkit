===================
gravity_disturbance
===================

- Calculates the gravity disturbance at a given latitude and longitude

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.gravity_disturbance import gravity_disturbance
    delta_g_h = gravity_disturbance(lat, lon, h, 'WGS84', clm, slm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/gravity_disturbance.py

.. autofunction:: geoid_toolkit.gravity_disturbance
