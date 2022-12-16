==============
norm_potential
==============

- Calculates the normal potential at a given latitude and height

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.norm_potential import norm_potential
    U, dU_dr, dU_dtheta = norm_potential(lat, lon, h, 'WGS84', lmax)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/norm_potential.py

.. autofunction:: geoid_toolkit.norm_potential

.. autofunction:: geoid_toolkit.norm_potential.cosine_even_zonals
