#!/usr/bin/env python
u"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (09/2020)
Reads the coefficients for a given gravity model file

GFZ International Centre for Global Earth Models (ICGEM)
    http://icgem.gfz-potsdam.de/

INPUTS:
    model_file: full path to *.gfc file with spherical harmonic coefficients

OPTIONS:
    LMAX: maximum degree and order of output spherical harmonic coefficients
    TIDE: tide system of output gravity fields
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential
    FLAG: string denoting data lines
    ZIP: input gravity field file is compressed in an archive file

OUTPUTS:
    l: spherical harmonic degree to maximum degree of model
    m: spherical harmonic order to maximum degree of model
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
    Updated 09/2021: define int/float precision to prevent deprecation warning
        update tidal offset to be able to change to and from any reference
        output spherical harmonic degree and order in dict
    Updated 03/2021: made degree of truncation LMAX a keyword argument
    Updated 07/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import io
import zipfile
import numpy as np
from geoid_toolkit.calculate_tidal_offset import calculate_tidal_offset

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(model_file, LMAX=None, TIDE=None,
    FLAG='gfc', ZIP=False):
    """
    Extract gravity model spherical harmonics from GFZ ICGEM gfc files

    Arguments
    ---------
    model_file: full path to gfc spherical harmonic data file
    LMAX: maximum degree and order of output spherical harmonics

    Keyword arguments
    -----------------
    TIDE: tide system of output gravity fields
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential
    FLAG: string denoting data lines
    ZIP: input gravity field file is compressed in an archive file

    Returns
    -------
    l: spherical harmonic degree to maximum degree of model
    m: spherical harmonic order to maximum degree of model
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

    #-- read data from compressed or gfc file
    if ZIP:
        #-- extract zip file with gfc file
        with zipfile.ZipFile(os.path.expanduser(model_file)) as zs:
            #-- find gfc file within zipfile
            gfc, = [io.BytesIO(zs.read(s)) for s in zs.namelist()
                if s.endswith('gfc')]
            #-- read input gfc data file
            file_contents = gfc.read().decode('ISO-8859-1').splitlines()
    else:
        #-- read input gfc data file
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
    LMAX = np.int64(model_input['max_degree']) if not LMAX else LMAX
    #-- output dimensions
    model_input['l'] = np.arange(LMAX+1)
    model_input['m'] = np.arange(LMAX+1)
    #-- allocate for each coefficient
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
        l1 = int(line_contents[1])
        m1 = int(line_contents[2])
        #-- if degree and order are below the truncation limits
        if ((l1 <= LMAX) and (m1 <= LMAX)):
            model_input['clm'][l1,m1] = np.float64(line_contents[3])
            model_input['slm'][l1,m1] = np.float64(line_contents[4])
            #-- check if model contains errors
            try:
                model_input['eclm'][l1,m1] = np.float64(line_contents[5])
                model_input['eslm'][l1,m1] = np.float64(line_contents[6])
            except:
                pass
    #-- calculate the tidal offset if changing the tide system
    if TIDE in ('mean_tide','zero_tide','tide_free'):
        #-- earth parameters
        GM = np.float64(model_input['earth_gravity_constant'])
        R = np.float64(model_input['radius'])
        model_input['clm'][2,0] += calculate_tidal_offset(TIDE,GM,R,'WGS84',
            REFERENCE=model_input['tide_system'])
        #-- update attribute for tide system
        model_input['tide_system'] = TIDE
    #-- return the spherical harmonics and parameters
    return model_input
