=========================
read_ICGEM_geoid_grids.py
=========================

- Reads geoid height spatial grids from the `GFZ International Centre for Global Earth Models (ICGEM) <http://icgem.gfz-potsdam.de/>`_
- Outputs spatial grids as netCDF4 files

Calling Sequence
################

.. code-block:: bash

    python read_ICGEM_geoid_grids.py --verbose FILE

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/scripts/read_ICGEM_geoid_grids.py

Inputs
######

1. ``FILE``: input geoid height spatial grids (``.gdf``)

Command Line Options
####################

- ``-F X``, ``--filename X``: Output netCDF4 filename (default: input with ``.nc`` suffix)
- ``-H X``, ``--header X``: Marker denoting the end of the header text
- ``-S X``, ``--spacing X``: Change output grid spacing
- ``-V``, ``--verbose``: Output information for each output file
- ``-M X``, ``--mode X``: Permission mode of output files
