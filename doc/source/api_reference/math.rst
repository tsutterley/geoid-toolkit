======================
``geoid_toolkit.math``
======================

- Special functions of mathematical physics

    * Computes the Gaussian weights as a function of degree using a normalized version of Christopher Jekeli's Gaussian averaging function
    * Computes fully normalized Legendre polynomials for an array of :math:`x` values and their first derivative
    * Creates pre-conditioned arrays for computing Clenshaw summations


Calling Sequence
================

.. code-block:: python

    import geoid_toolkit as geoidtk
    wl = 2.0*np.pi*geoidtk.math.gauss_weights(hw,LMAX)

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/math.py


General Methods
===============

.. autofunction:: geoid_toolkit.math.gauss_weights

.. autofunction:: geoid_toolkit.math.legendre_polynomials

.. autofunction:: geoid_toolkit.math.clenshaw_s_m

.. autofunction:: geoid_toolkit.math.clenshaw_ds_m

.. autofunction:: geoid_toolkit.math.clenshaw_ds_m_dr

.. autofunction:: geoid_toolkit.math.condon_shortley

.. autofunction:: geoid_toolkit.math.kronecker_delta
