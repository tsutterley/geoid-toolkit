==========================
corrected_geoid_undulation
==========================

- Calculates the topographically corrected geoidal undulation at a given latitude and longitude using an iterative approach

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.corrected_geoid_undulation import corrected_geoid_undulation
    N = corrected_geoid_undulation(lat, lon, 'WGS84', clm, slm, tclm, tslm, lmax, R, GM)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/corrected_geoid_undulation.py

.. autofunction:: geoid_toolkit.corrected_geoid_undulation
