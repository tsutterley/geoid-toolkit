#!/usr/bin/env python
"""
read_ICGEM_geoid_grids.py
Written by Tyler Sutterley (08/2026)
Reads geoid height spatial grids from the GFZ Geoid Calculation Service
    http://icgem.gfz-potsdam.de/home
Outputs spatial grids as netCDF4 files

INPUTS:
    input geoid height spatial grids (*.gdf)

COMMAND LINE OPTIONS:
    -F X, --filename X: Output filename (default: input with netCDF4 suffix)
    -H X, --header X: Marker denoting the end of the header text
    -S X, --spacing X: Change output grid spacing (via binning)
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/

UPDATE HISTORY:
    Updated 08/2026: use structured netCDF4 output to reduce redundancy
    Updated 06/2025: use import_dependency to import optional packages
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of geoid toolkit
        place netCDF4 import within try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: include utf-8 encoding in reads to be windows compliant
        check if gravity field data file is present in file-system
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: updated comments and argparse help text
    Updated 12/2020: using argparse to set parameters
    Updated 04/2019: verify that the divide count is greater than zero
    Updated 03/2018: generalized program using getopt to set parameters
        can regrid the data via binning with the --spacing option
        can output with a different filename with the --filename option
    Written 07/2014
"""

from __future__ import print_function

import sys
import re
import logging
import pathlib
import argparse
import datetime
import numpy as np
import geoid_toolkit as geoidtk


# PURPOSE: Reads .gdf grids from the GFZ calculation service
def read_ICGEM_geoid_grids(
    FILE,
    FILENAME=None,
    MARKER='',
    SPACING=None,
    VERBOSE=False,
    MODE=0o775,
):
    # create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # split filename into basename and extension
    FILE = pathlib.Path(FILE).expanduser().absolute()
    if not FILE.exists():
        raise FileNotFoundError(f'{str(FILE)} not found')
    # open input file and read contents
    with FILE.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()
    # number of lines contained in the file
    file_lines = len(file_contents)

    # counts the number of lines in the header
    count = 0
    # Reading over header text and extracting parameters
    HEADER = False
    parameters = {}
    while HEADER is False:
        # file line at count
        line = file_contents[count]
        if len(line) > 1:
            col = line.split()
            parameters[col[0].strip()] = col[1].strip()
        # find MARKER within line to set HEADER flag to True when found
        HEADER = bool(re.match(MARKER, line))
        # add 1 to counter
        count += 1

    # clean up dictionary of parameters
    for key in ['[deg.]', 'longitude', MARKER]:
        parameters = removekey(parameters, key)

    # extract necessary parameters
    latlimit_north = np.float64(parameters['latlimit_north'])
    latlimit_south = np.float64(parameters['latlimit_south'])
    longlimit_west = np.float64(parameters['longlimit_west'])
    longlimit_east = np.float64(parameters['longlimit_east'])
    # change grid spacing by binning data
    if SPACING is None:
        nlat = np.int64(parameters['latitude_parallels'])
        nlon = np.int64(parameters['longitude_parallels'])
        dlon = np.float64(parameters['gridstep'])
        dlat = np.float64(parameters['gridstep'])
    else:
        dlon, dlat = SPACING
        parameters['gridstep'] = f'{dlon:g},{dlat:g}'
        nlat = np.abs((latlimit_north - latlimit_south) / dlat).astype('i') + 1
        nlon = np.abs((longlimit_west - longlimit_east) / dlon).astype('i') + 1
        parameters['latitude_parallels'] = f'{nlat:d}'
        parameters['longitude_parallels'] = f'{nlon:d}'

    # output dataset
    dinput = {}
    # variable name and fill value
    functional = parameters.pop('functional')
    gapvalue = parameters.pop('gapvalue', np.nan)
    # allocate for output variable and mask
    dinput[functional] = np.ma.zeros((nlat, nlon), fill_value=gapvalue)
    dinput[functional].mask = np.zeros((nlat, nlon), dtype=bool)
    # create arrays of longitude and latitude
    dinput['lon'] = longlimit_west + np.arange(nlon) * dlon
    dinput['lat'] = latlimit_north - np.arange(nlat) * dlat
    # dictionary describing the output netCDF4 structure
    struct = dict(
        dimensions=('lat', 'lon'),
        variables={
            functional: ('lat', 'lon'),
        },
    )

    # get attributes from parameters
    attributes = {'ROOT': {}, functional: {}}
    # variable attributes
    variable_attributes = [
        'unit',
        'weighted_mean',
        'maxvalue',
        'minvalue',
        'signal_wrms',
        'zero_degree_term',
    ]
    for att_name in variable_attributes:
        attributes[functional][att_name] = parameters.pop(att_name, '')
    # attributes for longitude and latitude
    long_lat_unit = parameters.pop('unit_long_lat', 'degrees')
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
    # file-level attributes
    attributes['ROOT'].update(parameters)
    reference = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attributes['ROOT']['reference'] = reference

    # for each file line
    bin_count = np.zeros((nlat, nlon))
    for j in range(count, file_lines):
        col = np.array(file_contents[j].split(), dtype=np.float64)
        # calculating the lon/lat indice
        ilon = int((longlimit_west + col[0]) / dlon)
        ilat = int((latlimit_north - col[1]) / dlat)
        # if wanting data lat/lon
        dinput[functional][ilat, ilon] += np.float64(col[2])
        bin_count[ilat, ilon] += 1.0

    # take the mean of the binned data (if not regridding will divide by 1)
    ii, jj = np.nonzero(bin_count > 0)
    dinput[functional].data[ii, jj] /= bin_count[ii, jj]
    ii, jj = np.nonzero(bin_count == 0)
    dinput[functional].data[ii, jj] = dinput[functional].fill_value
    dinput[functional].mask[ii, jj] = True

    # output data and parameters to netCDF4
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    geoidtk.spatial.to_netCDF4(
        dinput,
        parameters,
        filename=FILENAME,
        structure=struct,
        data_type='structured',
    )
    # change permissions mode to MODE
    FILENAME.chmod(mode=MODE)


# PURPOSE: remove keys from an input dictionary
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads geoid height spatial grids from the ICGEM
            Geoid Calculation Service
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
    # marker denoting the end of the header text
    parser.add_argument(
        '--header',
        '-H',
        type=str,
        default='end_of_head',
        help='Marker denoting the end of the header text',
    )
    # change the output grid spacing by binning
    parser.add_argument(
        '--spacing',
        '-S',
        type=float,
        default=None,
        nargs=2,
        help='Change output grid spacing',
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
    read_ICGEM_geoid_grids(
        args.gravity,
        FILENAME=args.filename,
        MARKER=args.header,
        SPACING=args.spacing,
        VERBOSE=args.verbose,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
