=============================
compute_geoidal_undulation.py
=============================

- Computes geoid undulations from a gravity model for an input file (ascii, netCDF4, HDF5, geotiff)

Calling Sequence
################

.. code-block:: python

    python compute_geoidal_undulation.py --gravity <path_to_gravity_model> input_file output_file

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/scripts/compute_geoidal_undulation.py

Inputs
######

1. ``input_file``: name of input file
2. ``output_file``: name of output file

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
- ``--format X``: input and output data format

    * ``'csv'`` (default)
    * ``'netCDF4'``
    * ``'HDF5'``
    * ``'geotiff'``
- ``--variables X``: variable names of data in csv, HDF5 or netCDF4 file

    * for csv files: the order of the columns within the file
    * for HDF5 and netCDF4 files: time, y, x and data variable names
- ``-H X``, ``--header X``: number of header lines for csv files
- ``-t X``, ``--type X``: input data type

    * ``'drift'``: drift buoys or satellite/airborne altimetry (time per data point)
    * ``'grid'``: spatial grids or images (single time for all data points)
- ``--projection X``: spatial projection as EPSG code or PROJ4 string

    * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
- ``-V``, ``--verbose``: Verbose output of processing run
- ``-M X``, ``--mode X``: Permission mode of output file
