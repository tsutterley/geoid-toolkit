=====================
topographic_potential
=====================

- Calculates the potential at a given latitude and height using coefficients from a topographic model

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.topographic_potential import topographic_potential
    T = topographic_potential(lat, lon, refell, clm, slm, lmax, R, density)


`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/topographic_potential.py

.. autofunction:: geoid_toolkit.topographic_potential
