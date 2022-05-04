======================
compute_geoid_grids.py
======================

- Computes geoid undulations from a gravity model

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/scripts/compute_geoid_grids.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/compute_geoid_grids.py
    :func: arguments
    :prog: compute_geoid_grids.py
    :nodescription:
    :nodefault:

    --tide -T : @replace
        `Tide system of output geoid <http://mitgcm.org/~mlosch/geoidcookbook/node9.html>`_

        * ``'tide_free'``: no permanent direct and indirect tidal potentials
        * ``'mean_tide'``: permanent tidal potentials (direct and indirect)
        * ``'zero_tide'``: permanent direct tidal potential removed

    --projection -P : @after
        * ``4326``: latitude and longitude coordinates on WGS84 reference ellipsoid
