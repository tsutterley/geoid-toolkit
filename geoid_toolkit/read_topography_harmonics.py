#!/usr/bin/env python
u"""
read_topography_harmonics.py
Written by Tyler Sutterley (04/2022)
Reads the coefficients for a given topographic model file
http://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/

INPUTS:
    model_file: full path to file with spherical harmonic coefficients

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
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2017
"""
import numpy as np

#-- PURPOSE: read Earth 2014 topography harmonics
#-- http://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/
def read_topography_harmonics(model_file):
    """
    Reads `Earth 2014
    <https://ddfe.curtin.edu.au/gravitymodels/Earth2014/potential_model/readme_earth2014_potential_models.txt>`_
    topography harmonics from [Rexer2016]_

    Parameters
    ----------
    model_file: str
        full path to file with spherical harmonic coefficients

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
    density: float
        density of the Earth for the topography model

    References
    ----------
    .. [Rexer2016] M. Rexer, C. Hirt, S. Claessens, and R. Tenzer,
        "Layer-based modelling of the Earth's gravitational potential
        up to 10km scale in spherical harmonics in spherical and
        ellipsoidal approximation", *Surveys in Geophysics*, (2016).
        `doi:10.1007/s10712-016-9382-2 <https://doi.org/10.1007/s10712-016-9382-2>`_
    """
    dinput = np.fromfile(model_file, dtype=np.dtype('<f8'))
    #-- extract minimum and maximum spherical harmonic degree
    header = 2
    input_lmin,input_lmax = dinput[:header].astype(np.int)
    #-- number of spherical harmonic records for Clm and Slm
    n_down = ((input_lmin-1)**2 + 3*(input_lmin-1))/2 + 1
    n_up = (input_lmax**2 + 3*input_lmax)/2 + 1
    n_harm = n_up - n_down
    #-- dictionary of model parameters and output Ylms
    model_input = {}
    model_input['modelname'] = 'EARTH2014'
    model_input['density'] = 2670.0
    #-- extract cosine and sine harmonics
    ii,jj = np.tril_indices(input_lmax+1)
    #-- output dimensions
    model_input['l'] = np.arange(input_lmax+1)
    model_input['m'] = np.arange(input_lmax+1)
    model_input['clm'] = np.zeros((input_lmax+1,input_lmax+1))
    model_input['slm'] = np.zeros((input_lmax+1,input_lmax+1))
    model_input['clm'][ii,jj] = dinput[header:(header+n_harm)]
    model_input['slm'][ii,jj] = dinput[(header+n_harm):(header+2*n_harm)]
    return model_input