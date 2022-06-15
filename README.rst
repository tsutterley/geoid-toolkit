=============
geoid-toolkit
=============

|Language|
|License|
|Documentation Status|
|Binder|
|Pangeo|
|zenodo|

.. |Language| image:: https://img.shields.io/badge/python-v3.8-green.svg
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/github/license/tsutterley/geoid-toolkit
   :target: https://github.com/tsutterley/geoid-toolkit/blob/main/LICENSE

.. |PyPI Version| image:: https://img.shields.io/pypi/v/geoid-toolkit.svg
   :target: https://pypi.python.org/pypi/geoid-toolkit/

.. |Documentation Status| image:: https://readthedocs.org/projects/geoid-toolkit/badge/?version=latest
   :target: https://read-grace-harmonics.readthedocs.io/projects/geoid-toolkit/en/latest/?badge=latest

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/tsutterley/geoid-toolkit/main

.. |Pangeo| image:: https://img.shields.io/static/v1.svg?logo=Jupyter&label=PangeoBinderAWS&message=us-west-2&color=orange
   :target: https://aws-uswest2-binder.pangeo.io/v2/gh/tsutterley/geoid-toolkit/main?urlpath=lab

.. |zenodo| image:: https://zenodo.org/badge/199943310.svg
   :target: https://zenodo.org/badge/latestdoi/199943310

Python tools for obtaining and working with static gravity field coefficients
and calculating geoid heights

- `GFZ International Centre for Global Earth Models (ICGEM) <http://icgem.gfz-potsdam.de>`_

Dependencies
############

- `cartopy: Python package designed for geospatial data processing <https://scitools.org.uk/cartopy>`_
- `gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL) <https://pypi.python.org/pypi/GDAL>`_
- `h5py: Python interface for Hierarchal Data Format 5 (HDF5) <https://www.h5py.org/>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_
- `matplotlib: Python 2D plotting library <https://matplotlib.org>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/>`_
- `numpy: Scientific Computing Tools For Python <https://www.numpy.org>`_
- `pyproj: Python interface to PROJ library <https://pypi.org/project/pyproj/>`_

References
##########

    Drewes, Kuglitsch, Ad\ |aacute|\ m and R\ |oacute|\ zsa "The Geodesist's Handbook 2016",
    Journal of Geodesy, 90, 907-1205 (2016).
    `doi:10.1007/s00190-016-0948-z <https://doi.org/10.1007/s00190-016-0948-z>`_

    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005).
    `doi:10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_

    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
    the Recursive Computation of Very High Degree and Order Normalised
    Associated Legendre Functions", Journal of Geodesy (2002).
    `doi:10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_

    Ince, Barthelmes, Rei\ |szlig|\ land, Elger, F\ |ouml|\ rste, Flechtner, and Schuh,
    "ICGEM – 15 years of successful collection and distribution of global
    gravitational models, associated services, and future plans"
    Earth System Science Data, 11, 647–674 (2019).
    `doi:10.5194/essd-11-647-2019 <https://doi.org/10.5194/essd-11-647-2019>`_

    Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", (1981).

    Moazezi and Zomorrodian, "GGMCalc a software for calculation of the geoid
    undulation and the height anomaly using the iteration method, and
    classical gravity anomaly", Earth Science Informatics (2012).
    `doi:10.1007/s12145-012-0102-2 <https://doi.org/10.1007/s12145-012-0102-2>`_

    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
    Bollettino di Geodesia e Scienze, (1982)

    Wahr, Molenaar and Frank, "Time variability of the Earth's gravity field:
    Hydrological and oceanic effects and their possible detection using
    GRACE", Journal of Geophysical Research: Solid Earth, 103(B12),
    30205-30229, `doi:10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

Download
########

| The program homepage is:
| https://github.com/tsutterley/geoid-toolkit
| A zip archive of the latest version is available directly at:
| https://github.com/tsutterley/geoid-toolkit/archive/main.zip

Disclaimer
##########

This project contains work and contributions from the `scientific community <./CONTRIBUTORS.rst>`_.
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.
It is provided here for your convenience but *with no guarantees whatsoever*.

License
#######

The content of this project is licensed under the `Creative Commons Attribution 4.0 Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ and the source code is licensed under the `MIT license <LICENSE>`_.

.. |aacute|    unicode:: U+00E1 .. LATIN SMALL LETTER A WITH ACUTE
.. |oacute|    unicode:: U+00F3 .. LATIN SMALL LETTER O WITH ACUTE
.. |szlig|    unicode:: U+00DF .. LATIN SMALL LETTER SHARP S
.. |ouml|    unicode:: U+00F6 .. LATIN SMALL LETTER O WITH DIAERESIS
