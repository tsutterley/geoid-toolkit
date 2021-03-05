#!/usr/bin/env python
u"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (07/2020)
Reads the coefficients for a given gravity model file

GFZ International Centre for Global Earth Models (ICGEM)
    http://icgem.gfz-potsdam.de/

INPUTS:
    model_file: full path to *.gfc file with spherical harmonic coefficients

OPTIONS:
    LMAX: maximum degree and order of output spherical harmonic coefficients
    TIDE: tide system of output geoid
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
            this is the default (leaving the model as is)
        mean_tide: restores permanent tidal potentials (direct and indirect)
        zero_tide: restores permanent direct tidal potential
    FLAG: string denoting data lines

OUTPUTS:
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for the gravity model
    radius: semi-major axis of the Earth for the gravity model
    max_degree: maximum degree and order for the gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model (mean_tide, zero_tide, tide_free)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system

UPDATE HISTORY:
    Updated 03/2021: made degree of truncation LMAX a keyword argument
    Updated 07/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import numpy as np
from geoid_toolkit.calculate_tidal_offset import calculate_tidal_offset

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(model_file, LMAX=None, TIDE='tide_free', FLAG='gfc'):
    """
    Extract gravity model spherical harmonics from GFZ ICGEM gfc files

    Arguments
    ---------
    model_file: full path to *.gfc file with spherical harmonic coefficients
    LMAX: maximum degree and order of output spherical harmonic coefficients

    Keyword arguments
    -----------------
    TIDE: tide system of output geoid
    FLAG: string denoting data lines

    Returns
    -------
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for gravity model
    radius: semi-major axis of the Earth for gravity model
    max_degree: maximum degree and order for gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model
    """

    #-- read input data
    with open(os.path.expanduser(model_file),'r') as f:
        file_contents = f.read().splitlines()
    #-- python dictionary with model input and headers
    model_input = {}
    #-- extract parameters from header
    header_parameters = ['modelname','earth_gravity_constant','radius',
        'max_degree','errors','norm','tide_system']
    parameters_regex = '(' + '|'.join(header_parameters) + ')'
    header = [l for l in file_contents if re.match(parameters_regex,l)]
    for line in header:
        #-- split the line into individual components
        line_contents = line.split()
        model_input[line_contents[0]] = line_contents[1]
    #-- set degree of truncation from model if not presently set
    LMAX = np.int(model_input['max_degree']) if not LMAX else LMAX
    #-- allocate for each Coefficient
    model_input['clm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['slm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['eclm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['eslm'] = np.zeros((LMAX+1,LMAX+1))
    #-- reduce file_contents to input data using data marker flag
    input_data = [l for l in file_contents if re.match(FLAG,l)]
    #-- for each line of data in the gravity file
    for line in input_data:
        #-- split the line into individual components replacing fortran d
        line_contents = re.sub('d','e',line,flags=re.IGNORECASE).split()
        #-- degree and order for the line
        l1 = np.int(line_contents[1])
        m1 = np.int(line_contents[2])
        #-- if degree and order are below the truncation limits
        if ((l1 <= LMAX) and (m1 <= LMAX)):
            model_input['clm'][l1,m1] = np.float(line_contents[3])
            model_input['slm'][l1,m1] = np.float(line_contents[4])
            model_input['eclm'][l1,m1] = np.float(line_contents[5])
            model_input['eslm'][l1,m1] = np.float(line_contents[6])
    #-- calculate the tidal offset if changing the tide system
    if TIDE in ('mean_tide','zero_tide'):
        model_input['tide_system'] = TIDE
        GM = np.float(model_input['earth_gravity_constant'])
        R = np.float(model_input['radius'])
        model_input['clm'][2,0] += calculate_tidal_offset(TIDE,GM,R,'WGS84')
    #-- return the spherical harmonics and parameters
    return model_input
