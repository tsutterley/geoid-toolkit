============
Installation
============

``geoid-toolkit`` is available for download from the `GitHub repository <https://github.com/tsutterley/geoid-toolkit>`_,
the `Python Package Index (pypi) <https://pypi.org/project/geoid-toolkit/>`_,
and from `conda-forge <https://anaconda.org/conda-forge/geoid-toolkit>`_.


The simplest installation for most users will likely be using ``conda`` or ``mamba``:

.. code-block:: bash

    conda install -c conda-forge geoid-toolkit

``conda`` installed versions of ``geoid-toolkit`` can be upgraded to the latest stable release:

.. code-block:: bash

    conda update geoid-toolkit

To use the development repository, please fork ``geoid-toolkit`` into your own account and then clone onto your system:

.. code-block:: bash

    git clone https://github.com/tsutterley/geoid-toolkit.git

``geoid-toolkit`` can then be installed within the package directory using ``pip``:

.. code-block:: bash

    python3 -m pip install --user .

The development version of ``geoid-toolkit`` can also be installed directly from GitHub using ``pip``:

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/geoid-toolkit.git
