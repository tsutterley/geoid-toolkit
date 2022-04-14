======================
compute_geoid_grids.py
======================

- Computes geoid undulations from a gravity model

Calling Sequence
################

.. code-block:: python

    python compute_geoid_grids.py --gravity <path_to_gravity_model> output_file

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/scripts/compute_geoid_grids.py

Inputs
######

1. ``output_file``: name of output file

Command Line Options
####################

- ``-D X``, ``--directory X``: Working data directory
- ``-G X``, ``--gravity X``: Gravity model file to use (.gfc format)
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree (level of truncation)
- ``-T X``, ``--tide X``: `tide system of output geoid <http://mitgcm.org/~mlosch/geoidcookbook/node9.html>`_

    * ``'tide_free'``: no permanent direct and indirect tidal potentials
    * ``'mean_tide'``: permanent tidal potentials (direct and indirect)
    * ``'zero_tide'``: permanent direct tidal potential removed
- ``-R X``, ``--radius X``: Gaussian smoothing radius (km)
- ``--format X``: output data format

    * ``'csv'`` (default)
    * ``'netCDF4'``
    * ``'HDF5'``
    * ``'geotiff'``
- ``-S X``, ``--spacing X``: Output grid spacing
- ``-B X``, ``--bounds X``: output grid extents [xmin,xmax,ymin,ymax]
- ``--projection X``: spatial projection as EPSG code or PROJ4 string

    * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
- ``-V``, ``--verbose``: Verbose output of processing run
- ``-M X``, ``--mode X``: Permission mode of output file
