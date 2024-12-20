=============================
compute_geoidal_undulation.py
=============================

- Computes geoid undulations from a gravity model for an input file (ascii, netCDF4, HDF5, GTiff, cog)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/scripts/compute_geoidal_undulation.py

Calling Sequence
################

.. argparse::
    :filename: compute_geoidal_undulation.py
    :func: arguments
    :prog: compute_geoidal_undulation.py
    :nodescription:
    :nodefault:

    --tide -T : @replace
        `tide system of output geoid <http://mitgcm.org/~mlosch/geoidcookbook/node9.html>`_

        * ``'tide_free'``: no permanent direct and indirect tidal potentials
        * ``'mean_tide'``: permanent tidal potentials (direct and indirect)
        * ``'zero_tide'``: permanent direct tidal potential removed

    --variables : @after
        * for csv files: the order of the columns within the file
        * for HDF5 and netCDF4 files: time, y, x and data variable names

    --type -t : @after
        * ``'drift'``: drift buoys or satellite/airborne altimetry (time per data point)
        * ``'grid'``: spatial grids or images (single time for all data points)

    --projection -P : @after
        * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
