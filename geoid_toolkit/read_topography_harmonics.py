#!/usr/bin/env python
"""
read_topography_harmonics.py
Written by Tyler Sutterley (08/2026)
Reads the coefficients for a given topographic model file
http://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/

INPUTS:
    model_file: full path to file with spherical harmonic coefficients
    LMAX: maximum degree and order of output spherical harmonics

OUTPUTS:
    l: spherical harmonic degree to maximum degree of model
    m: spherical harmonic order to maximum degree of model
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the topography model
    density: density of the Earth for the topography model

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 08/2026: add LMAX option to truncate/extend spherical harmonics
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2017
"""

import pathlib
import numpy as np


__all__ = [
    'read_topography_harmonics',
]


# PURPOSE: read Earth 2014 topography harmonics
# http://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/
def read_topography_harmonics(model_file: str | pathlib.Path, **kwargs):
    """
    Reads Earth 2014 topography harmonics from :cite:t:`Rexer:2016gr`

    Parameters
    ----------
    model_file: str
        full path to file with spherical harmonic coefficients
    LMAX: int or NoneType, default None
        maximum degree and order of output spherical harmonics

    Returns
    -------
    l: int
        spherical harmonic degree of model
    m: int
        spherical harmonic order to maximum degree of model
    clm: float
        cosine spherical harmonics of topographic data
    slm: float
        sine spherical harmonics of topographic data
    modelname: str
        name of the topography model
    max_degree: str
        Maximum degree and order for gravity model
    density: float
        density of the Earth for the topography model
    """
    dinput = np.fromfile(model_file, dtype=np.dtype('<f8'))
    # extract minimum and maximum spherical harmonic degree
    header = 2
    input_lmin, input_lmax = dinput[:header].astype(np.int64)
    # number of spherical harmonic records for Clm and Slm
    n_down = ((input_lmin - 1) ** 2 + 3 * (input_lmin - 1)) / 2 + 1
    n_up = (input_lmax**2 + 3 * input_lmax) / 2 + 1
    n_harm = n_up - n_down
    # set degree of truncation from model if not presently set
    LMAX = kwargs.get('LMAX', None) or input_lmax
    # dictionary of model parameters and output Ylms
    model_input = {}
    # set maximum degree attribute
    model_input['max_degree'] = str(LMAX)
    # set model name and density attribute
    model_input['modelname'] = 'EARTH2014'
    model_input['density'] = 2670.0
    # extract cosine and sine harmonics
    ii, jj = np.tril_indices(input_lmax + 1)
    # output dimensions
    model_input['l'] = np.arange(input_lmax + 1)
    model_input['m'] = np.arange(input_lmax + 1)
    model_input['clm'] = np.zeros((input_lmax + 1, input_lmax + 1))
    model_input['slm'] = np.zeros((input_lmax + 1, input_lmax + 1))
    model_input['clm'][ii, jj] = dinput[header : (header + n_harm)]
    model_input['slm'][ii, jj] = dinput[
        (header + n_harm) : (header + 2 * n_harm)
    ]
    # truncate or extend spherical harmonics to maximum degree
    # if LMAX is equal to the input lmax: then use original harmonics
    if LMAX < input_lmax:
        # truncate spherical harmonics to maximum degree
        model_input['l'] = model_input['l'][: LMAX + 1]
        model_input['m'] = model_input['m'][: LMAX + 1]
        model_input['clm'] = model_input['clm'][: LMAX + 1, : LMAX + 1]
        model_input['slm'] = model_input['slm'][: LMAX + 1, : LMAX + 1]
    elif LMAX > input_lmax:
        # extend spherical harmonics to maximum degree
        clm = np.zeros((LMAX + 1, LMAX + 1))
        slm = np.zeros((LMAX + 1, LMAX + 1))
        clm[: input_lmax + 1, : input_lmax + 1] = model_input['clm']
        slm[: input_lmax + 1, : input_lmax + 1] = model_input['slm']
        # update model_input dictionary with extended spherical harmonics
        model_input['l'] = np.arange(LMAX + 1)
        model_input['m'] = np.arange(LMAX + 1)
        model_input['clm'] = clm
        model_input['slm'] = slm
    # return spherical harmonics and parameters
    return model_input
