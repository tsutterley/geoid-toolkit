============================
read_topography_harmonics.py
============================

- Reads the coefficients for a given `topographic model file <http://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/>`_

Calling Sequence
################

.. code-block:: python

    from geoid_toolkit.read_topography_harmonics import read_topography_harmonics
    Ylms = read_topography_harmonics(model_file)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/read_topography_harmonics.py

Arguments
#########

1. ``model_file``: full path to file with spherical harmonic coefficients

Returns
#######

- ``l``: spherical harmonic degree to maximum degree of model
- ``m``: spherical harmonic order to maximum degree of model
- ``clm``: cosine spherical harmonics of input data
- ``slm``: sine spherical harmonics of input data
- ``eclm``: cosine spherical harmonic standard deviations of type errors
- ``eslm``: sine spherical harmonic standard deviations of type errors
- ``modelname``: name of the topography model
- ``density``: density of the Earth for the topography model
