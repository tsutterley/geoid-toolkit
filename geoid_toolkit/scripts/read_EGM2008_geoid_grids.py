#!/usr/bin/env python
"""
read_EGM2008_geoid_grids.py
Written by Tyler Sutterley (08/2026)
Reads EGM2008 geoid height spatial grids from unformatted binary files
    provided by the National Geospatial-Intelligence Agency
Outputs spatial grids as netCDF4 files

NGA Office of Geomatics
    https://earth-info.nga.mil/

INPUTS:
    input 2.5x2.5 arcminute geoid height spatial grids

COMMAND LINE OPTIONS:
    -F X, --filename X: Output filename
    -n X, --love X: Degree 2 load Love number
    -G, --gzip: Input file is gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/

UPDATE HISTORY:
    Updated 08/2026: use structured netCDF4 output to reduce redundancy
    Updated 07/2026: input files can be gzip compressed
        added output attribute for the degree 2 load Love number
    Updated 06/2025: use import_dependency to import optional packages
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of geoid toolkit
        place netCDF4 import within try/except statements
    Written 06/2022
"""

from __future__ import print_function

import sys
import gzip
import logging
import pathlib
import argparse
import datetime
import numpy as np
import geoid_toolkit as geoidtk


def read_EGM2008_geoid_grids(
    FILE,
    FILENAME=None,
    LOVE=0.3,
    GZIP=False,
    VERBOSE=False,
    MODE=0o775,
):
    # create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # check that data file is present in file system
    FILE = pathlib.Path(FILE).expanduser().absolute()
    if not FILE.exists():
        raise FileNotFoundError(f'{str(FILE)} not found')
    opener = gzip.open if GZIP else open
    with opener(FILE, mode='rb') as f:
        # open input file and read contents
        file_contents = np.frombuffer(f.read(), dtype='<f4')

    # set grid parameters
    dlon, dlat = (2.5 / 60.0), (2.5 / 60.0)
    latlimit_north, latlimit_south = (90.0, -90.0)
    longlimit_west, longlimit_east = (0.0, 360.0)
    # boundary parameters
    nlat = np.abs((latlimit_north - latlimit_south) / dlat).astype('i') + 1
    nlon = np.abs((longlimit_west - longlimit_east) / dlon).astype('i') + 1

    # variable and file-level attributes
    attributes = dict(ROOT={})
    fill_value = -9999.0
    # root attributes
    attributes['ROOT']['source'] = 'EGM2008'
    attributes['ROOT']['reference'] = 'http://earth-info.nima.mil/GandG/'
    attributes['ROOT']['earth_gravity_constant'] = 0.3986004415e15
    attributes['ROOT']['radius'] = 0.63781363e07
    attributes['ROOT']['max_degree'] = 2190
    attributes['ROOT']['norm'] = 'fully_normalized'
    reference = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['reference'] = reference
    # latitude and longitude
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lon']['valid_min'] = longlimit_west
    attributes['lon']['valid_max'] = longlimit_east
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    attributes['lat']['valid_min'] = latlimit_south
    attributes['lat']['valid_max'] = latlimit_north
    # geoid_h
    attributes['geoid_h'] = {}
    attributes['geoid_h']['long_name'] = 'Geoidal_Undulation'
    attributes['geoid_h']['description'] = (
        'Geoid undulations with respect to WGS84'
    )
    attributes['geoid_h']['units'] = 'meters'
    attributes['geoid_h']['fill_value'] = fill_value
    attributes['geoid_h']['tide_system'] = 'tide_free'
    attributes['geoid_h']['source'] = 'EGM2008'
    # geoid_free2mean
    attributes['geoid_free2mean'] = {}
    attributes['geoid_free2mean']['long_name'] = 'Free-to-Mean conversion'
    attributes['geoid_free2mean']['description'] = (
        'Additive value to convert geoid heights from the tide-free '
        'system to the mean-tide system'
    )
    attributes['geoid_free2mean']['units'] = 'meters'
    attributes['geoid_free2mean']['fill_value'] = fill_value
    attributes['geoid_free2mean']['tide_system'] = 'tide_free'
    attributes['geoid_free2mean']['source'] = 'derived'
    attributes['geoid_free2mean']['k2'] = LOVE

    # dictionary describing the output netCDF4 structure
    struct = dict(
        dimensions=('lat', 'lon'),
        variables={
            'geoid_h': ('lat', 'lon'),
            'geoid_free2mean': ('lat', 'lon'),
        },
    )

    # output variables
    dinput = {}
    # create arrays of longitude and latitude
    dinput['lon'] = longlimit_west + np.arange(nlon) * dlon
    dinput['lat'] = latlimit_north - np.arange(nlat) * dlat
    # geoid undulation
    dinput['geoid_h'] = np.ma.zeros((nlat, nlon), dtype=np.float32)
    dinput['geoid_h'].fill_value = fill_value
    # reshape data to matrix
    geoid_h = file_contents.reshape(nlat, nlon + 1)
    dinput['geoid_h'][:, :-1] = geoid_h[:, 1:-1]
    # repeat values for 360
    dinput['geoid_h'][:, -1] = dinput['geoid_h'][:, 0]

    # calculate Legendre polynomial of degree 2 (unnormalized)
    gridlon, gridlat = np.meshgrid(dinput['lon'], dinput['lat'])
    P2 = 0.5 * (3.0 * np.sin(np.radians(gridlat)) ** 2 - 1.0)
    # offset for converting from tide_free to mean_tide
    # from Rapp 1991 (Consideration of Permanent Tidal Deformation)
    dinput['geoid_free2mean'] = np.ma.zeros((nlat, nlon), dtype=np.float32)
    dinput['geoid_free2mean'].fill_value = fill_value
    dinput['geoid_free2mean'][:, :] = -0.198 * P2 * (1.0 + LOVE)

    # output data and parameters to netCDF4
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    geoidtk.spatial.to_netCDF4(
        dinput,
        attributes,
        filename=FILENAME,
        structure=struct,
        data_type='structured',
    )
    # change permissions mode to MODE
    FILENAME.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads EGM2008 geoid height spatial grids from
            unformatted binary files provided by the National
            Geospatial-Intelligence Agency
            """
    )
    # command line parameters
    parser.add_argument(
        'gravity',
        type=pathlib.Path,
        help='Geoid height spatial grid file',
    )
    # output filename (will default to input file with netCDF4 suffix)
    parser.add_argument(
        '--filename',
        '-F',
        type=pathlib.Path,
        help='Output netCDF4 filename',
    )
    # load love number of degree 2 (default EGM2008 value)
    parser.add_argument(
        '--love',
        '-n',
        type=float,
        default=0.3,
        help='Degree 2 load Love number',
    )
    # input file is gzip compressed
    parser.add_argument(
        '--gzip',
        '-G',
        default=False,
        action='store_true',
        help='Input file is gzip compressed',
    )
    # verbose will output information about each output file
    parser.add_argument(
        '--verbose',
        '-V',
        default=False,
        action='store_true',
        help='Output information for each output file',
    )
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permissions mode of output files',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # verify input and output files
    args.gravity = pathlib.Path(args.gravity).expanduser().absolute()
    # set output file from input filename if not entered
    if not args.filename:
        args.filename = args.gravity.with_name(f'{args.gravity.stem}.nc')
    else:
        args.filename = pathlib.Path(args.filename).expanduser().absolute()

    # run program
    read_EGM2008_geoid_grids(
        args.gravity,
        FILENAME=args.filename,
        LOVE=args.love,
        GZIP=args.gzip,
        VERBOSE=args.verbose,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
