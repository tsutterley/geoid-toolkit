#!/usr/bin/env python
"""
math.py
Written by Tyler Sutterley (08/2026)
Special functions of mathematical physics

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    Hofmann-Wellenhof and Moritz (2005)
        "Physical Geodesy"
        http://www.springerlink.com/content/978-3-211-33544-4
    Christopher Jekeli (1981)
        "Alternative Methods to Smooth the Earth's Gravity Field"
        http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf

NOTES:
    IDL code gauss_weights.pro was written by Sean Swenson

    Differences from recurs function in combine.mac.f:
    weighting from gauss_weights is normalized outside of the function
        wt = 2.0*pi*gauss_weights(rad,LMAX)
    weighting from recurs is normalized inside of the function
    call recurs(alpha,bcoef) calculates bcoef up to LMAX 150 (=wt[0:150])
        alpha = alog(2.)/(1.-cos(rad/6371.))

UPDATE HISTORY:
    Updated 08/2026: separated mathematical functions into a new math module
        simplify Clenshaw summation to reduce memory usage
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use complex form of spherical harmonics for summations
        use np.radians to convert from degrees to radians
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: added option for setting minimum value threshold
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 09/2020: verify dimensions of x variable
    Updated 08/2020: prevent zero divisions by changing u==0 to eps of data type
    Updated 07/2020: added function docstrings
    Updated 07/2017: added Gaussian smoothing with option GAUSS
        changed dtypes to long double for high degree and order models
        added first derivative of Legendre polynomials (dpl)
        added option ASTYPE to output as different variable types
    Updated 06/2015: adjusted threshold from 1e-9 to 1e-10
    Updated 12/2014: updated comments and header text updating full reference
    Updated 02/2014: changed variables from ints to floats to prevent truncation
    Written 03/2013
"""

from __future__ import annotations
import numpy as np


__all__ = [
    'gauss_weights',
    'legendre_polynomials',
    'clenshaw_s_m',
    'clenshaw_ds_m',
    'clenshaw_ds_m_dr',
    'condon_shortley',
    'kronecker_delta',
]


def gauss_weights(hw, LMAX, CUTOFF=1e-10):
    r"""
    Computes the Gaussian weights as a function of degree using
    a normalized form from :cite:t:`Jekeli:1981vj`

    Parameters
    ----------
    hw: float
        Gaussian smoothing radius in kilometers

        Radius :math:`r` corresponds to the distance at which the
        weight drops to half its peak value at the shortest wavelength
    LMAX: int
        Maximum degree of spherical harmonic coefficients
    CUTOFF: float, default 1e-10
        minimum value for tail of Gaussian averaging function

    Returns
    -------
    wl: float
        degree dependent weighting function
    """
    # allocate for output weights
    wl = np.zeros((LMAX + 1))
    # radius of the Earth in km
    rad_e = 6371.0
    if hw < CUTOFF:
        # distance is smaller than cutoff
        wl[:] = 1.0 / (2.0 * np.pi)
    else:
        # calculate gaussian weights using recursion
        b = np.log(2.0) / (1.0 - np.cos(hw / rad_e))
        # weight for degree 0
        wl[0] = 1.0 / (2.0 * np.pi)
        # weight for degree 1
        wl[1] = wl[0] * (
            (1.0 + np.exp(-2.0 * b)) / (1.0 - np.exp(-2.0 * b)) - 1.0 / b
        )
        # valid flag
        valid = True
        # spherical harmonic degree
        l = 2
        # while valid (within cutoff)
        # and spherical harmonic degree is less than LMAX
        while valid and (l <= LMAX):
            # calculate weight with recursion
            wl[l] = (1.0 - 2.0 * l) / b * wl[l - 1] + wl[l - 2]
            # weight is less than cutoff
            if wl[l] < CUTOFF:
                # set all weights to cutoff
                wl[l : LMAX + 1] = CUTOFF
                # set valid flag
                valid = False
            # add 1 to l
            l += 1
    # return the gaussian weights
    return wl


def legendre_polynomials(lmax, x, ASTYPE=np.float64):
    r"""
    Computes fully-normalized Legendre polynomials and their first derivative
    following :cite:t:`HofmannWellenhof:2006hy`

    Calculates Legendre polynomials for zonal harmonics (order 0)

    Parameters
    ----------
    lmax: int
        maximum degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically :math:`\cos(\theta)`, where :math:`\theta`
        is the colatitude in radians
    ASTYPE: np.dtype, default np.float64
        output variable data type

    Returns
    -------
    pl: np.ndarray
        fully-normalized Legendre polynomials
    dpl: np.ndarray
        first derivative of Legendre polynomials
    """
    # verify dimensions
    x = np.atleast_1d(x).flatten().astype(ASTYPE)
    # size of the x array
    nx = len(x)
    # verify data type of spherical harmonic truncation
    lmax = np.int64(lmax)
    # output matrix of normalized legendre polynomials
    pl = np.zeros((lmax + 1, nx), dtype=ASTYPE)
    # output matrix of First derivative of Legendre polynomials
    dpl = np.zeros((lmax + 1, nx), dtype=ASTYPE)
    # dummy matrix for the recurrence relation
    ptemp = np.zeros((lmax + 1, nx), dtype=ASTYPE)

    # u is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    # for x=cos(th): u=sin(th)
    u = np.sqrt(1.0 - x**2)
    # update where u==0 to eps of data type to prevent invalid divisions
    u0 = np.flatnonzero(u == 0)
    u[u0] = np.finfo(u.dtype).eps

    # Initialize the recurrence relation
    # ptemp is a dummy array of length lmax+1 storing unnormalized values
    ptemp[0, :] = 1.0
    ptemp[1, :] = x
    # Normalization is geodesy convention
    pl[0, :] = ptemp[0, :]
    pl[1, :] = np.sqrt(3.0) * ptemp[1, :]
    for l in range(2, lmax + 1):
        ptemp[l, :] = (((2.0 * l) - 1.0) / l) * x * ptemp[l - 1, :] - (
            (l - 1.0) / l
        ) * ptemp[l - 2, :]
        # Normalization is geodesy convention
        pl[l, :] = np.sqrt((2.0 * l) + 1.0) * ptemp[l, :]
        # Overwrite polar case (x == +/-1)
        pl[l, u0] = np.sqrt((2.0 * l) + 1.0) * x[u0] ** l

    # First derivative of Legendre polynomials
    for l in range(1, lmax + 1):
        fl = np.sqrt(((l**2.0) * (2.0 * l + 1.0)) / (2.0 * l - 1.0))
        dpl[l, :] = (1.0 / u) * (l * x * pl[l, :] - fl * pl[l - 1, :])

    # return the legendre polynomials and their first derivative
    return (pl, dpl)


# PURPOSE: compute Clenshaw summation of the fully normalized associated
# Legendre's function for constant order m
def clenshaw_s_m(
    t: np.ndarray,
    q: np.ndarray,
    m: int,
    Ylm1: np.ndarray,
    lmax: int,
    SCALE: float = 1e-280,
):
    r"""
    Compute conditioned arrays for Clenshaw summation from the fully-normalized
    associated Legendre's function for an order m

    Parameters
    ----------
    t: np.ndarray
        :math:`\cos(\theta)`, where :math:`\theta` is the colatitude in radians
    q: np.ndarray
        degree dependent factors to apply
    m: int
        spherical harmonic order
    Ylm1: np.ndarray
        complex form of spherical harmonics
    lmax: int
        maximum spherical harmonic degree (truncation limit)
    SCALE: float, default 1e-280
        scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    cs_m: np.ndarray
        conditioned array for clenshaw summation
    """
    # allocate for output matrix
    N = len(t)
    cs_m = np.zeros((N), dtype=np.clongdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.clongdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        cs_m[:] = np.copy(ylm[lmax, lmax])
    elif m == (lmax - 1):
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        cs_m[:] = a_lm * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:] = np.copy(s_mm_l)
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = (
            t * q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1) * (ll + 1))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        cs_m[:] = np.copy(s_mm_l)
    # return rescaled cs_m
    return cs_m / SCALE


# PURPOSE: compute Clenshaw summation of derivative with respect to latitude
# of the fully normalized associated Legendre's function for constant order m
def clenshaw_ds_m(t, u, q, m, Ylm1, lmax, SCALE=1e-280):
    r"""
    Compute Clenshaw summation of derivative with respect to latitude of
    the fully normalized associated Legendre's function for order m

    Parameters
    ----------
    t: np.ndarray
        :math:`\cos(\theta)`, where :math:`\theta` is the colatitude in radians
    u: np.ndarray
        :math:`\sin(\theta)`, where :math:`\theta` is the colatitude in radians
    q: np.ndarray
        degree dependent factors to apply
    m: int
        spherical harmonic order
    Ylm1: np.ndarray
        complex form of spherical harmonics
    lmax: int
        maximum spherical harmonic degree (truncation limit)
    SCALE: float, default 1e-280
        scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    dcs_m: np.ndarray
        conditioned array for clenshaw summation
    """
    # allocate for output matrix
    N = len(t)
    dcs_m = np.zeros((N), dtype=np.clongdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.clongdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        dcs_m[:] = mm * t * u * ylm[lmax, lmax]
    elif m == (lmax - 1):
        a_lm = q * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        s_dot_mm = a_lm * ylm[lmax, lmax - 1]
        s_mm_l = a_lm * t * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
        dcs_m[:] = mm * t * u * s_mm_l - u * s_dot_mm
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = q * np.sqrt(
            ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
        )
        s_dot_mm_minus_2 = 0.0
        s_dot_mm_minus_1 = a_lm * s_mm_minus_2
        s_mm_minus_1 = a_lm * t * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = q * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_dot_mm = (
                a_lm * (s_dot_mm_minus_1 * t + s_mm_minus_1)
                - b_lm * s_dot_mm_minus_2
            )
            s_mm_l = a_lm * t * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
            s_dot_mm_minus_2 = np.copy(s_dot_mm_minus_1)
            s_dot_mm_minus_1 = np.copy(s_dot_mm)
        dcs_m[:] = mm * t * u * s_mm_l - u * s_dot_mm
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        s_dot_mm_minus_1 = a_lm * s_mm_minus_2
        s_mm_minus_1 = a_lm * t * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = q * np.sqrt(
                ((2.0 * ll + 1.0) * (2.0 * ll + 3.0)) / ((ll + 1) * (ll + 1))
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            s_dot_mm = (
                a_lm * (s_dot_mm_minus_1 * t + s_mm_minus_1)
                - b_lm * s_dot_mm_minus_2
            )
            s_mm_l = a_lm * t * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
            s_dot_mm_minus_2 = np.copy(s_dot_mm_minus_1)
            s_dot_mm_minus_1 = np.copy(s_dot_mm)
        dcs_m[:] = -u * s_dot_mm
    # return rescaled dcs_m
    return dcs_m / SCALE


# PURPOSE: compute Clenshaw summation of derivative with respect to radius of
# the fully normalized associated Legendre's function for constant order m
def clenshaw_ds_m_dr(t, q, m, Ylm1, lmax, SCALE=1e-280):
    r"""
    Compute Clenshaw summation of derivative with respect to radius of
    the fully normalized associated Legendre's function for order m

    Parameters
    ----------
    t: np.ndarray
        :math:`\cos(\theta)`, where :math:`\theta` is the colatitude in radians
    q: np.ndarray
        degree dependent factors to apply
    m: int
        spherical harmonic order
    Ylm1: np.ndarray
        complex form of spherical harmonics
    lmax: int
        maximum spherical harmonic degree (truncation limit)
    SCALE: float, default 1e-280
        scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    dcs_m_dr: np.ndarray
        conditioned array for clenshaw summation
    """
    # allocate for output matrix
    N = len(t)
    dcs_m_dr = np.zeros((N), dtype=np.clongdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.clongdouble)
    # convert lmax and m to float
    lm = np.longdouble(lmax)
    mm = np.longdouble(m)
    if m == lmax:
        dcs_m_dr[:] = -(lm + 1.0) * (ylm[lmax, lmax])
    elif m == (lmax - 1):
        ds_mm_dr_c_pre_1 = -(lm + 1.0) * ylm[lmax, lmax - 1]
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        dcs_m_dr[:] = (
            a_lm * ds_mm_dr_c_pre_1 * ylm[lmax, lmax - 1]
            - lm * ylm[lmax - 1, lmax - 1]
        )
    elif (m <= (lmax - 2)) and (m >= 1):
        ds_mm_dr_minus_2 = -(lm + 1.0) * ylm[lmax, m]
        a_lm = (
            t
            * q
            * np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
        )
        ds_mm_dr_minus_1 = a_lm * ds_mm_dr_minus_2 - lm * ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            ds_mm_dr = (
                a_lm * ds_mm_dr_minus_1
                - b_lm * ds_mm_dr_minus_2
                - (ll + 1.0) * ylm[l, m]
            )
            ds_mm_dr_minus_2 = np.copy(ds_mm_dr_minus_1)
            ds_mm_dr_minus_1 = np.copy(ds_mm_dr)
        dcs_m_dr[:] = np.copy(ds_mm_dr)
    elif m == 0:
        ds_mm_dr_minus_2 = -(lm + 1.0) * ylm[lmax, 0]
        a_lm = (
            t * q * np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm))
        )
        ds_mm_dr_minus_1 = a_lm * ds_mm_dr_minus_2 - lm * ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.longdouble(l)
            a_lm = (
                t
                * q
                * np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1) * (ll + 1))
                )
            )
            b_lm = np.power(q, 2) * np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2) * (ll + 2) * (2.0 * ll + 1.0))
            )
            ds_mm_dr = (
                a_lm * ds_mm_dr_minus_1
                - b_lm * ds_mm_dr_minus_2
                - (ll + 1.0) * ylm[l, 0]
            )
            ds_mm_dr_minus_2 = np.copy(ds_mm_dr_minus_1)
            ds_mm_dr_minus_1 = np.copy(ds_mm_dr)
        dcs_m_dr[:] = np.copy(ds_mm_dr)
    # return rescaled dcs_m_dr
    return dcs_m_dr / SCALE


def condon_shortley(m: int | np.ndarray):
    r"""
    Computes the Condon-Shortley phase :math:`(-1)^m` for order :math:`m`

    Parameters
    ----------
    m: int or np.ndarray
        Order of the Legendre polynomials
    """
    return np.power(-1.0, m)


def kronecker_delta(
    i: int | np.ndarray,
    j: int | np.ndarray,
):
    r"""
    Computes the Kronecker delta :math:`\delta_{ij}` function

    .. math::
        \delta_{ij} =
            \begin{cases}
                1 & \text{if } i = j \\
                0 & \text{if } i \neq j
            \end{cases}

    Parameters
    ----------
    i: int or np.ndarray
        First index
    j: int or np.ndarray
        Second index
    """
    return 1.0 * (i == j)
