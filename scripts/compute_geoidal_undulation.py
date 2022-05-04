#!/usr/bin/env python
u"""
compute_geoidal_undulation.py
Written by Tyler Sutterley (05/2022)
Computes geoid undulations from a gravity model for an input file

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates
    geotiff file with bands in spatial coordinates

COMMAND LINE OPTIONS:
    -G X, --gravity X: Gravity model file to use (.gfc format)
    -l X, --lmax X: maximum spherical harmonic degree (level of truncation)
    -T X, --tide X: tide system of output geoid
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
    -R X, --radius X: Gaussian smoothing radius (km)
    --format X: input and output data format
        csv (default)
        netCDF4
        HDF5
        geotiff
    --variables X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5 and netCDF4 files: time, y, x and data variable names
    -H X, --header X: number of header lines for csv files
    -t X, --type X: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
    -P X, --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permission mode of output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    spatial: utilities for reading, writing and operating on spatial data
    utilities.py: download and management utilities for syncing files
    geoid_undulation.py: geoidal undulation at a given latitude and longitude
    read_ICGEM_harmonics.py: reads the coefficients for a given gravity model file
    real_potential.py: real potential at a latitude and height for gravity model
    norm_potential.py: normal potential of an ellipsoid at a latitude and height
    norm_gravity.py: normal gravity of an ellipsoid at a latitude and height
    ref_ellipsoid.py: Computes parameters for a reference ellipsoid
    gauss_weights.py: Computes Gaussian weights as a function of degree

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 11/2021: add function for attempting to extract projection
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 11/2020: options to read from and write to geotiff image files
    Updated 10/2020: using argparse to set command line parameters
    Updated 09/2020: can use HDF5 and netCDF4 as inputs and outputs
    Written 07/2017
"""
from __future__ import print_function

import sys
import os
import pyproj
import logging
import argparse
import numpy as np
import geoid_toolkit.spatial
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from geoid_toolkit.geoid_undulation import geoid_undulation
from geoid_toolkit.utilities import convert_arg_line_to_args

#-- PURPOSE: try to get the projection information for the input file
def get_projection(attributes, PROJECTION):
    #-- coordinate reference system string from file
    try:
        crs = pyproj.CRS.from_string(attributes['projection'])
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- EPSG projection code
    try:
        crs = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- no projection can be made
    raise pyproj.exceptions.CRSError

#-- PURPOSE: read csv, netCDF or HDF5 data and compute geoid undulation
def compute_geoidal_undulation(model_file, input_file, output_file,
    LMAX=None,
    TIDE=None,
    GAUSS=0,
    FORMAT='csv',
    VARIABLES=[],
    HEADER=0,
    TYPE='drift',
    PROJECTION='4326',
    VERBOSE=False,
    MODE=0o775):

    #-- create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- read gravity model Ylms and change tide if specified
    Ylms = read_ICGEM_harmonics(model_file, LMAX=LMAX, TIDE=TIDE)
    R = np.float64(Ylms['radius'])
    GM = np.float64(Ylms['earth_gravity_constant'])
    LMAX = np.int64(Ylms['max_degree'])
    #-- reference to WGS84 ellipsoid
    REFERENCE = 'WGS84'
    #-- invalid value
    fill_value = -9999.0
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- file level attributes
    attrib['ROOT']['model'] = Ylms['modelname']
    attrib['ROOT']['ellipsoid'] = REFERENCE
    attrib['ROOT']['earth_gravity_constant'] = GM
    attrib['ROOT']['radius'] = R
    #-- latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude'
    attrib['lat']['units'] = 'Degrees_North'
    #-- longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude'
    attrib['lon']['units'] = 'Degrees_East'
    #-- geoid undulation
    attrib['geoid_h'] = {}
    attrib['geoid_h']['units'] = 'm'
    attrib['geoid_h']['long_name'] = 'Geoidal_Undulation'
    args = (Ylms['modelname'],Ylms['max_degree'])
    attrib['geoid_h']['description'] = ('{0}_geoidal_undulation_'
        'computed_from_degree_{1}_gravity_model.').format(*args)
    attrib['geoid_h']['tide_system'] = Ylms['tide_system']
    attrib['geoid_h']['degree_of_truncation'] = LMAX
    attrib['geoid_h']['_FillValue'] = fill_value

    #-- read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        dinput = geoid_toolkit.spatial.from_ascii(input_file,
            columns=VARIABLES,
            header=HEADER)
    elif (FORMAT == 'netCDF4'):
        dinput = geoid_toolkit.spatial.from_netCDF4(input_file,
            timename=VARIABLES[0],
            xname=VARIABLES[2],
            yname=VARIABLES[1],
            varname=VARIABLES[3])
    elif (FORMAT == 'HDF5'):
        dinput = geoid_toolkit.spatial.from_HDF5(input_file,
            timename=VARIABLES[0],
            xname=VARIABLES[2],
            yname=VARIABLES[1],
            varname=VARIABLES[3])
    elif (FORMAT == 'geotiff'):
        dinput = geoid_toolkit.spatial.from_geotiff(input_file)
        #-- copy global geotiff attributes for projection and grid parameters
        for att_name in ['projection','wkt','spacing','extent']:
            attrib[att_name] = dinput['attributes'][att_name]

    #-- converting x,y from projection to latitude/longitude
    crs1 = get_projection(dinput['attributes'], PROJECTION)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (TYPE == 'grid'):
        ny,nx = (len(dinput['y']),len(dinput['x']))
        gridx,gridy = np.meshgrid(dinput['x'],dinput['y'])
        lon,lat = transformer.transform(gridx,gridy)
        #-- calculate geoid at coordinates and reshape to output
        N = geoid_undulation(lat.flatten(), lon.flatten(),
            REFERENCE, Ylms['clm'], Ylms['slm'],
            LMAX, R, GM, GAUSS=GAUSS).reshape(ny,nx)
    elif (TYPE == 'drift'):
        lon,lat = transformer.transform(dinput['x'],dinput['y'])
        #-- calculate geoid at coordinates
        N = geoid_undulation(lat, lon,
            REFERENCE, Ylms['clm'], Ylms['slm'],
            LMAX, R, GM, GAUSS=GAUSS)
    #-- replace fill values with fill value
    N = np.ma.asarray(N,dtype=np.float64)
    N.mask = np.copy(dinput['data'].mask)
    N.fill_value = np.copy(fill_value)
    N.data[N.mask] = N.fill_value

    #-- output to file
    output = {'lon':lon,'lat':lat,'geoid_h':N}
    if (FORMAT == 'csv'):
        geoid_toolkit.spatial.to_ascii(output, attrib, output_file,
            delimiter=',', columns=['lat','lon','geoid_h'])
    elif (FORMAT == 'netCDF4'):
        geoid_toolkit.spatial.to_netCDF4(output, attrib, output_file)
    elif (FORMAT == 'HDF5'):
        geoid_toolkit.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        geoid_toolkit.spatial.to_geotiff(output, attrib, output_file,
            varname='geoid_h')
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates geoidal undulations for an input file
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = convert_arg_line_to_args
    #-- command line options
    #-- input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input file to run')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Computed output file')
    #-- set gravity model file to use
    parser.add_argument('--gravity','-G',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Gravity model file to use')
    #-- maximum spherical harmonic degree (level of truncation)
    parser.add_argument('--lmax','-l',
        type=int, help='Maximum spherical harmonic degree')
    #-- tide system of output geoid
    #-- tide_free: no permanent direct and indirect tidal potentials
    #-- mean_tide: permanent tidal potentials (direct and indirect)
    #-- zero_tide: permanent direct tidal potential removed
    parser.add_argument('--tide','-T',
        type=str, default='tide_free',
        choices=['tide_free','mean_tide','zero_tide'],
        help='Tide system of output geoid')
    #-- Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    #-- input and output data format
    parser.add_argument('--format','-F',
        type=str, default='csv', choices=('csv','netCDF4','HDF5','geotiff'),
        help='Input and output data format')
    #-- variable names (for csv names of columns)
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['time','lat','lon','data'],
        help='Variable names of data in input file')
    #-- number of header lines for csv files
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines for csv files')
    #-- input data type
    #-- drift: drift buoys or satellite/airborne altimetry (time per data point)
    #-- grid: spatial grids or images (single time for all data points)
    parser.add_argument('--type','-t',
        type=str, default='drift',
        choices=('drift','grid'),
        help='Input data type')
    #-- spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    #-- verbose output of processing run
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of output file')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- set output file from input filename if not entered
    if not args.outfile:
        model,_ = os.path.splitext(os.path.basename(args.gravity))
        fileBasename,fileExtension = os.path.splitext(args.infile)
        vars = (fileBasename,model,fileExtension)
        args.outfile = '{0}_{1}{2}'.format(*vars)

    #-- run geoid undulation program for input file
    compute_geoidal_undulation(args.gravity, args.infile, args.outfile,
        LMAX=args.lmax,
        TIDE=args.tide,
        GAUSS=args.radius,
        FORMAT=args.format,
        VARIABLES=args.variables,
        HEADER=args.header,
        TYPE=args.type,
        PROJECTION=args.projection,
        VERBOSE=args.verbose,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
