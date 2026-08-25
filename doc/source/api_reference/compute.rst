=========================
``geoid_toolkit.compute``
=========================

- Utilities for computing functionals from a gravity model

    * Geoidal undulation at a given latitude and longitude using an iterative approach
    * Gravity anomaly at a given latitude and longitude using different methods
    * Gravity disturbance at a given latitude and longitude
    * Height anomaly at a given latitude and longitude using an iterative approach
    * Real potential at a given latitude and height using coefficients from a gravity model
    * Potential at a given latitude and height using coefficients from a topographic model


Calling Sequence
================

.. code-block:: python

    import geoid_toolkit as geoidtk
    N = geoidtk.compute.geoid_height(lon, lat, gravity_model_file, ELLIPSOID='WGS84')

`Source code`__

.. __: https://github.com/polargeodesy/geoid-toolkit/blob/main/geoid_toolkit/compute.py


General Methods
===============

.. autofunction:: geoid_toolkit.compute.geoid_height

.. autofunction:: geoid_toolkit.compute.geoid_undulation

.. autofunction:: geoid_toolkit.compute.corrected_geoid_undulation

.. autofunction:: geoid_toolkit.compute.gravity_anomaly

.. autofunction:: geoid_toolkit.compute.gravity_disturbance

.. autofunction:: geoid_toolkit.compute.height_anomaly

.. autofunction:: geoid_toolkit.compute.real_potential

.. autofunction:: geoid_toolkit.compute.topographic_potential
