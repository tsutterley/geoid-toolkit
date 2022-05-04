#!/usr/bin/env python
u"""
compute_geoid_grids.py
Written by Tyler Sutterley (05/2022)
Computes geoid undulations from a gravity model

COMMAND LINE OPTIONS:
    -G X, --gravity X: Gravity model file to use (.gfc format)
    -l X, --lmax X: maximum spherical harmonic degree (level of truncation)
    -T X, --tide X: tide system of output geoid
        http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
    -R X, --radius X: Gaussian smoothing radius (km)
    --format X: output data format
        csv (default)
        netCDF4
        HDF5
        geotiff
    -S X, --spacing X: output grid spacing
    -B X, --bounds X: output grid extents [xmin,xmax,ymin,ymax]
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
    Written 03/2022
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

#-- PURPOSE: try to get the projection information
def get_projection(PROJECTION):
    #-- EPSG projection code
    try:
        crs = pyproj.CRS.from_epsg(int(PROJECTION))
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

#-- PURPOSE: compute geoid undulation for a grid
def compute_geoid_grids(model_file, output_file,
    LMAX=None,
    TIDE=None,
    GAUSS=0,
    FORMAT='csv',
    BOUNDS=[],
    SPACING=[],
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

    #-- converting x,y from projection to latitude/longitude
    crs1 = get_projection(PROJECTION)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- dictionary of coordinate reference system variables
    crs_to_dict = crs1.to_dict()
    #-- Climate and Forecast (CF) Metadata Conventions
    if (crs1.to_epsg() == 4326):
        y_cf,x_cf = crs1.cs_to_cf()
    else:
        x_cf,y_cf = crs1.cs_to_cf()

    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    output,attrib = ({},{})
    #-- file level attributes
    attrib['ROOT'] = {}
    attrib['ROOT']['model'] = Ylms['modelname']
    attrib['ROOT']['ellipsoid'] = REFERENCE
    attrib['ROOT']['earth_gravity_constant'] = GM
    attrib['ROOT']['radius'] = R
    #-- projection attributes
    attrib['crs'] = {}
    #-- add projection attributes
    attrib['crs']['standard_name'] = crs1.to_cf()['grid_mapping_name'].title()
    attrib['crs']['spatial_epsg'] = crs1.to_epsg()
    attrib['crs']['spatial_ref'] = crs1.to_wkt()
    attrib['crs']['proj4_params'] = crs1.to_proj4()
    for att_name,att_val in crs1.to_cf().items():
        attrib['crs'][att_name] = att_val
    if ('lat_0' in crs_to_dict.keys() and (crs1.to_epsg() != 4326)):
        attrib['crs']['latitude_of_projection_origin'] = \
            crs_to_dict['lat_0']
    #-- x and y
    attrib['x'],attrib['y'] = ({},{})
    for att_name in ['long_name','standard_name','units']:
        attrib['x'][att_name] = x_cf[att_name]
        attrib['y'][att_name] = y_cf[att_name]
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

    #-- projection variable
    output['crs'] = np.array((),dtype=np.byte)
    #-- spacing and bounds of output grid
    dx,dy = np.broadcast_to(np.atleast_1d(SPACING),(2,))
    xmin,xmax,ymin,ymax = np.copy(BOUNDS)
    # create x and y from spacing and bounds
    output['x'] = np.arange(xmin+dx/2.0,xmax+dx,dx)
    if (FORMAT == 'geotiff'):
        output['y'] = np.arange(ymax-dx/2.0,ymin-dy,-dy)
    else:
        output['y'] = np.arange(ymin+dx/2.0,ymax+dy,dy)
    ny,nx = (len(output['y']),len(output['x']))
    gridx,gridy = np.meshgrid(output['x'],output['y'])
    lon,lat = transformer.transform(gridx,gridy)
    #-- calculate geoid at coordinates and reshape to output
    N = geoid_undulation(lat.flatten(), lon.flatten(),
        REFERENCE, Ylms['clm'], Ylms['slm'],
        LMAX, R, GM, GAUSS=GAUSS).reshape(ny,nx)

    #-- replace fill values with fill value
    output['geoid_h'] = np.ma.asarray(N,dtype=np.float64)
    output['geoid_h'].mask = np.logical_not(np.isfinite(N))
    output['geoid_h'].fill_value = np.copy(fill_value)
    output['geoid_h'].data[output['geoid_h'].mask] = fill_value

    #-- output to file
    if (FORMAT == 'csv'):
        geoid_toolkit.spatial.to_ascii(output, attrib, output_file,
            delimiter=',', columns=['y','x','geoid_h'])
    elif (FORMAT == 'netCDF4'):
        geoid_toolkit.spatial.to_netCDF4(output, attrib, output_file)
    elif (FORMAT == 'HDF5'):
        geoid_toolkit.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        #-- copy global geotiff attributes for projection and grid parameters
        attrib['wkt'] = crs1.to_wkt()
        attrib['spacing'] = (dx, -dy)
        attrib['extent'] = np.copy(BOUNDS)
        geoid_toolkit.spatial.to_geotiff(output, attrib, output_file,
            varname='geoid_h')
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Computes geoid undulations from a gravity model
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = convert_arg_line_to_args
    #-- command line options
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Output file')
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
    #-- Output data format
    parser.add_argument('--format','-F',
        type=str, default='csv', choices=('csv','netCDF4','HDF5','geotiff'),
        help='Output data format')
    #-- output grid spacing
    parser.add_argument('--spacing','-S',
        type=float, default=1.0, nargs='+',
        help='Output grid spacing')
    #-- bounds of output mosaic
    parser.add_argument('--bounds','-B', type=float,
        nargs=4, default=[-180.0,180.0,-90.0,90.0],
        metavar=('xmin','xmax','ymin','ymax'),
        help='Output grid extents')
    #-- spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    #-- verbose output of processing run
    #-- print information about each output file
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

    #-- run geoid grid program
    compute_geoid_grids(args.gravity, args.outfile,
        LMAX=args.lmax,
        TIDE=args.tide,
        GAUSS=args.radius,
        FORMAT=args.format,
        BOUNDS=args.bounds,
        SPACING=args.spacing,
        PROJECTION=args.projection,
        VERBOSE=args.verbose,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
