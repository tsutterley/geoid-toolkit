=======================
``geoid_toolkit.datum``
=======================

- Computes parameters for a reference ellipsoid
- Calculates the normal gravity of an ellipsoid at a given latitude and height
- Calculates the normal potential at a given latitude and height


Calling Sequence
================

.. code-block:: python

    import geoid_toolkit as geoidtk
    params = geoidtk.datum.ref_ellipsoid(refell)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/datum.py


General Methods
===============

.. autofunction:: geoid_toolkit.datum.ref_ellipsoid

.. autofunction:: geoid_toolkit.datum.norm_gravity

.. autofunction:: geoid_toolkit.datum.norm_potential
