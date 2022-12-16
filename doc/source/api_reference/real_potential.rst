==============
real_potential
==============

- Calculates the real potential at a given latitude and height using coefficients from a gravity model

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.real_potential import real_potential
    W, dW_dr, dW_dtheta = real_potential(lat, lon, h, clm, slm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/real_potential.py

.. autofunction:: geoid_toolkit.real_potential
