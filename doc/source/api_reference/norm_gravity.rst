===============
norm_gravity.py
===============

- Calculates the normal gravity of an ellipsoid at a given latitude and height and calculates the derivative with respect to height

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.norm_gravity import norm_gravity
    gamma_h,dgamma_dh = norm_gravity(latitude, height, 'WGS84')

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/norm_gravity.py

.. autofunction:: geoid_toolkit.norm_gravity
