================
``interpolate``
================

- Routines to interpolate data from pre-computed spatial grids

Calling Sequence
################

.. code-block:: python

    import geoid_toolkit.interpolate
    geoid_h = geoid_toolkit.interpolate.geoid_height(lon, lat, model='EGM2008')

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/interpolate.py

.. autoclass:: geoid_toolkit.interpolate.Interpolate
    :members:
    :private-members:

.. autofunction:: geoid_toolkit.interpolate.geoid_height

