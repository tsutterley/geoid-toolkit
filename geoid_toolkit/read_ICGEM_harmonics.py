#!/usr/bin/env python
"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (08/2026)
Reads the coefficients for a given gravity model file

GFZ International Centre for Global Earth Models (ICGEM)
    http://icgem.gfz-potsdam.de/

INPUTS:
    model_file: full path to *.gfc file with spherical harmonic coefficients

OPTIONS:
    LMAX: maximum degree and order of output spherical harmonic coefficients
    ELLIPSOID: reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967
        GRS80 = Geodetic Reference System 1980
        WGS72 = World Geodetic System 1972
        WGS84 = World Geodetic System 1984
        ATS77 = Quasi-earth centred ellipsoid for ATS77
        NAD27 = North American Datum 1927
        NAD83 = North American Datum 1983
        INTER = International
        KRASS = Krassovsky (USSR)
        MAIRY = Modified Airy (Ireland 1965/1975)
        TOPEX = TOPEX/POSEIDON ellipsoid
        EGM96 = EGM 1996 gravity model
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data
    TIDE: tide system of output gravity fields
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
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

UPDATE HISTORY:
    Updated 08/2026: move tidal offset function back into this program
        rename output dictionary to gfc for consistency with other functions
        raise ValueError for invalid tide system inputs
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        check if gravity field coefficients file is present in file-system
    Updated 10/2021: ellipsoid option for semi-major axis when changing tides
    Updated 09/2021: define int/float precision to prevent deprecation warning
        update tidal offset to be able to change to and from any reference
        output spherical harmonic degree and order in dict
    Updated 03/2021: made degree of truncation LMAX a keyword argument
    Updated 07/2020: added function docstrings
    Updated 07/2019: split read and wrapper function into separate files
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""

import re
import io
import pathlib
import zipfile
import numpy as np

__all__ = [
    'read_ICGEM_harmonics',
    '_tidal_offset',
]


# PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(
    model_file: str | pathlib.Path,
    **kwargs,
):
    """
    Extract gravity model spherical harmonics from GFZ ICGEM ``gfc`` files

    Parameters
    ----------
    model_file: str
        full path to gfc spherical harmonic data file
    LMAX: int or NoneType, default None
        maximum degree and order of output spherical harmonics
    ELLIPSOID: str
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model

    TIDE: str or NoneType, default None
        Permanent tide system of output gravity fields

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    FLAG: str, default 'gfc'
        Flag denoting data lines
    ZIP: bool, default False
        Gravity field file is compressed in an archive file

    Returns
    -------
    l: int
        spherical harmonic degree of model
    m: int
        spherical harmonic order of model
    clm: float
        cosine spherical harmonics of input data
    slm: float
        sine spherical harmonics of input data
    eclm: float
        cosine spherical harmonic standard deviations of type errors
    eslm: float
        sine spherical harmonic standard deviations of type errors
    modelname: str
        Name of the gravity model
    earth_gravity_constant: str
        GM constant of the Earth for gravity model
    radius: str
        Semi-major axis of the Earth for gravity model
    max_degree: str
        Maximum degree and order for gravity model
    errors: str
        Error type of the gravity model
    norm: str
        Normalization of the spherical harmonics
    tide_system: str
        Permanent tide system of gravity model
    """
    # set default keyword arguments
    kwargs.setdefault('ELLIPSOID', 'WGS84')
    kwargs.setdefault('TIDE', None)
    kwargs.setdefault('FLAG', 'gfc')
    kwargs.setdefault('ZIP', False)

    # tilde-expansion of input file
    if not isinstance(model_file, io.IOBase):
        model_file = pathlib.Path(model_file).expanduser().absolute()
        # check that data file is present in file system
        if not model_file.exists():
            raise FileNotFoundError(f'{str(model_file)} not found')

    # read data from compressed or gfc file
    if kwargs['ZIP']:
        # extract zip file with gfc file
        with zipfile.ZipFile(model_file) as zs:
            # find gfc file within zipfile
            (s,) = [s for s in zs.namelist() if s.endswith('gfc')]
            # read gravity field coefficients from zip
            f = io.BytesIO(zs.read(s))
            file_contents = f.read().decode('ISO-8859-1').splitlines()
    else:
        # read gravity field coefficients file
        with open(model_file, mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()

    # python dictionary with model input and headers
    gfc = {}
    # extract parameters from header
    header_parameters = [
        'modelname',
        'earth_gravity_constant',
        'radius',
        'max_degree',
        'errors',
        'norm',
        'tide_system',
    ]
    parameters_regex = '(' + '|'.join(header_parameters) + ')'
    header = [l for l in file_contents if re.match(parameters_regex, l)]
    for line in header:
        # split the line into individual components
        line_contents = line.split()
        gfc[line_contents[0]] = line_contents[1]
    # set degree of truncation from model if not presently set
    LMAX = kwargs.get('LMAX') or np.int64(gfc['max_degree'])
    # update maximum degree attribute if truncating
    if LMAX != np.int64(gfc['max_degree']):
        gfc['max_degree'] = str(LMAX)
    # output dimensions
    gfc['l'] = np.arange(LMAX + 1)
    gfc['m'] = np.arange(LMAX + 1)
    # allocate for each coefficient
    gfc['clm'] = np.zeros((LMAX + 1, LMAX + 1))
    gfc['slm'] = np.zeros((LMAX + 1, LMAX + 1))
    gfc['eclm'] = np.zeros((LMAX + 1, LMAX + 1))
    gfc['eslm'] = np.zeros((LMAX + 1, LMAX + 1))
    # reduce file_contents to input data using data marker flag
    input_data = [l for l in file_contents if re.match(kwargs['FLAG'], l)]
    # for each line of data in the gravity file
    for line in input_data:
        # split the line into individual components replacing fortran d
        line_contents = re.sub('d', 'e', line, flags=re.IGNORECASE).split()
        # degree and order for the line
        l1 = int(line_contents[1])
        m1 = int(line_contents[2])
        # if degree and order are below the truncation limits
        if (l1 <= LMAX) and (m1 <= LMAX):
            gfc['clm'][l1, m1] = np.float64(line_contents[3])
            gfc['slm'][l1, m1] = np.float64(line_contents[4])
            # check if model contains errors
            try:
                gfc['eclm'][l1, m1] = np.float64(line_contents[5])
                gfc['eslm'][l1, m1] = np.float64(line_contents[6])
            except Exception as exc:
                pass
    # calculate the tidal offset if changing the tide system
    if kwargs['TIDE'] in ('mean_tide', 'zero_tide', 'tide_free'):
        # earth parameters
        GM = np.float64(gfc['earth_gravity_constant'])
        R = np.float64(gfc['radius'])
        gfc['clm'][2, 0] += _tidal_offset(
            kwargs['TIDE'],
            GM,
            R,
            kwargs['ELLIPSOID'],
            REFERENCE=gfc['tide_system'],
        )
        # update attribute for tide system
        gfc['tide_system'] = kwargs['TIDE']
    # return the spherical harmonics and parameters
    return gfc


def _tidal_offset(
    TIDE: str,
    GM: float,
    R: float,
    refell: str,
    LOVE: float = 0.3,
    REFERENCE: str = 'tide_free',
):
    """
    Calculates the spherical harmonic offset to change permanent tide systems
    :cite:p:`HofmannWellenhof:2006hy,Losch:2003ve`

    Parameters
    ----------
    TIDE: str
        Output permanent tidal system

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    R: float
        Average radius used in gravity model
    GM: float
        Geocentric gravitational constant used in gravity model
    refell: str
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model
    LOVE: float, default 0.3
        Load love number for degree 2
    REFERENCE: str, default 'tide_free'
        Original permanent tidal system of gravity modeld

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed

    Returns
    -------
    delta: float
        Offset for changing to tide system
    """
    # import reference ellipsoid tools
    from geoid_toolkit.datum import ref_ellipsoid

    # get ellipsoid parameters for refell
    ellip = ref_ellipsoid(refell)
    # standard gravitational acceleration
    gamma = 9.80665
    trans = (-0.198 * gamma * R**3) / (np.sqrt(5.0) * GM * ellip['a'] ** 2)
    # conversion to switch to tide free
    if REFERENCE == 'tide_free':
        tide_free_conv = 0.0
    elif REFERENCE == 'mean_tide':
        tide_free_conv = -(1.0 + LOVE)
    elif REFERENCE == 'zero_tide':
        tide_free_conv = -LOVE
    else:
        raise ValueError(f'Invalid input permanent tide system {REFERENCE}')
    # conversion for each tidal system
    if TIDE == 'mean_tide':
        conv = (1.0 + LOVE) + tide_free_conv
    elif TIDE == 'zero_tide':
        conv = LOVE + tide_free_conv
    elif TIDE == 'tide_free':
        conv = 0.0 + tide_free_conv
    else:
        raise ValueError(f'Invalid output permanent tide system {TIDE}')
    # return the C20 offset to change tide systems
    delta = conv * trans
    return delta
