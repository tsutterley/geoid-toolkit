=======
spatial
=======

Utilities for reading, writing and operating on spatial data

 - Can read ascii, netCDF4, HDF5 or geotiff files
 - Can output to ascii, netCDF4, HDF5 or geotiff files

Calling Sequence
----------------

Reading a netCDF4 file

.. code-block:: python

    import geoid_toolkit.spatial
    dinput = geoid_toolkit.spatial.from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    import geoid_toolkit.spatial
    dinput = geoid_toolkit.spatial.from_HDF5(path_to_HDF5_file)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/spatial.py

General Methods
===============


.. autofunction:: geoid_toolkit.spatial.case_insensitive_filename

.. autofunction:: geoid_toolkit.spatial.data_type

.. autofunction:: geoid_toolkit.spatial.from_file

.. autofunction:: geoid_toolkit.spatial.from_ascii

.. autofunction:: geoid_toolkit.spatial.from_netCDF4

.. autofunction:: geoid_toolkit.spatial.from_HDF5

.. autofunction:: geoid_toolkit.spatial.from_geotiff

.. autofunction:: geoid_toolkit.spatial.to_ascii

.. autofunction:: geoid_toolkit.spatial.to_netCDF4

.. autofunction:: geoid_toolkit.spatial.to_HDF5

.. autofunction:: geoid_toolkit.spatial.to_geotiff

.. autofunction:: geoid_toolkit.spatial.expand_dims

.. autofunction:: geoid_toolkit.spatial.convert_ellipsoid

.. autofunction:: geoid_toolkit.spatial.compute_delta_h

.. autofunction:: geoid_toolkit.spatial.wrap_longitudes

.. autofunction:: geoid_toolkit.spatial.to_cartesian

.. autofunction:: geoid_toolkit.spatial.to_sphere

.. autofunction:: geoid_toolkit.spatial.to_geodetic

.. autofunction:: geoid_toolkit.spatial._moritz_iterative

.. autofunction:: geoid_toolkit.spatial._bowring_iterative

.. autofunction:: geoid_toolkit.spatial._zhu_closed_form

.. autofunction:: geoid_toolkit.spatial.geocentric_latitude
