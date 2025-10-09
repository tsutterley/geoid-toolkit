===========================
geoid-toolkit Documentation
===========================

Welcome to the documentation for ``geoid-toolkit``, a set of Python tools
for obtaining and working with static gravity field coefficients
and calculating geoid heights.

Introduction
------------

.. grid:: 2 2 2 2
    :padding: 0

    .. grid-item-card::  Installation
      :text-align: center
      :link: ./getting_started/Install.html

      :material-outlined:`download;5em`

    .. grid-item-card::  Examples
      :text-align: center
      :link: ./user_guide/Examples.html

      :material-outlined:`apps;5em`

Contribute
----------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Guidelines
      :text-align: center
      :link: ./getting_started/Contributing.html

      :material-outlined:`groups;5em`

    .. grid-item-card::  Code of Conduct
      :text-align: center
      :link: ./getting_started/Code-of-Conduct.html

      :material-outlined:`gavel;5em`

    .. grid-item-card::  Issues
      :text-align: center
      :link: https://github.com/tsutterley/geoid-toolkit/issues

      :material-outlined:`bug_report;5em`

    .. grid-item-card::  Citation Information
      :text-align: center
      :link: ./project/Citations.html

      :material-outlined:`alternate_email;5em`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    getting_started/Install.rst
    getting_started/Contributing.rst
    getting_started/Code-of-Conduct.rst
    getting_started/Resources.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: User Guide

    user_guide/Examples.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: API Reference

    api_reference/calculate_geoid_undulation.rst
    api_reference/calculate_tidal_offset.rst
    api_reference/corrected_geoid_undulation.rst
    api_reference/gauss_weights.rst
    api_reference/geoid_undulation.rst
    api_reference/gravity_anomaly.rst
    api_reference/gravity_disturbance.rst
    api_reference/height_anomaly.rst
    api_reference/legendre_polynomials.rst
    api_reference/norm_gravity.rst
    api_reference/norm_potential.rst
    api_reference/read_ICGEM_harmonics.rst
    api_reference/read_topography_harmonics.rst
    api_reference/real_potential.rst
    api_reference/ref_ellipsoid.rst
    api_reference/spatial.rst
    api_reference/topographic_potential.rst
    api_reference/utilities.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Use Cases

    api_reference/compute_geoid_grids.rst
    api_reference/compute_geoidal_undulation.rst
    api_reference/read_EGM2008_geoid_grids.rst
    api_reference/read_ICGEM_geoid_grids.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Project Details

    project/Contributors.rst
    project/Licenses.rst
    project/Testing.rst
    project/Citations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Bibliography

    project/Bibliography.rst
