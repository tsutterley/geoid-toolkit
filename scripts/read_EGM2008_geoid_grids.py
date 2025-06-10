#!/usr/bin/env python
u"""
read_EGM2008_geoid_grids.py
Written by Tyler Sutterley (06/2025)
Reads EGM2008 geoid height spatial grids from unformatted binary files
provided by the National Geospatial-Intelligence Agency
Outputs spatial grids as netCDF4 files

NGA Office of Geomatics
    https://earth-info.nga.mil/

INPUTS:
    input 2.5x2.5 arcminute geoid height spatial grids

COMMAND LINE OPTIONS:
    -F X, --filename X: Output filename
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/

UPDATE HISTORY:
    Updated 06/2025: use import_dependency to import optional packages
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of geoid toolkit
        place netCDF4 import within try/except statements
    Written 06/2022
"""
from __future__ import print_function

import logging
import pathlib
import argparse
import datetime
import numpy as np
import geoid_toolkit as geoidtk

# attempt imports
netCDF4 = geoidtk.utilities.import_dependency('netCDF4')

def read_EGM2008_geoid_grids(FILE, FILENAME=None, LOVE=0.3,
    VERBOSE=False, MODE=0o775):

    # create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # check that data file is present in file system
    FILE = pathlib.Path(FILE).expanduser().absolute()
    if not FILE.exists():
        raise FileNotFoundError(f'{str(FILE)} not found')
    # open input file and read contents
    file_contents = np.fromfile(FILE, dtype='<f4')

    # set grid parameters
    dlon,dlat = (2.5/60.0), (2.5/60.0)
    latlimit_north, latlimit_south = (90.0, -90.0)
    longlimit_west, longlimit_east = (0.0, 360.0)
    # boundary parameters
    nlat = np.abs((latlimit_north-latlimit_south)/dlat).astype('i') + 1
    nlon = np.abs((longlimit_west-longlimit_east)/dlon).astype('i') + 1

    # variable and file-level attributes
    attributes = dict(geoid_h={}, geoid_free2mean={}, ROOT={})
    # root attributes
    attributes['ROOT']['source'] = 'EGM2008'
    attributes['ROOT']['reference'] = 'http://earth-info.nima.mil/GandG/'
    attributes['ROOT']['earth_gravity_constant'] = 0.3986004415E+15
    attributes['ROOT']['radius'] = 0.63781363E+07
    attributes['ROOT']['max_degree'] = 2190
    attributes['ROOT']['norm'] = 'fully_normalized'
    # geoid_h
    attributes['geoid_h']['long_name'] = 'Geoidal_Undulation'
    attributes['geoid_h']['description'] = ('Geoid undulations with '
        'respect to WGS84')
    attributes['geoid_h']['units'] = 'meters'
    attributes['geoid_h']['fill_value'] = -9999.0
    attributes['geoid_h']['tide_system'] = 'tide_free'
    attributes['geoid_h']['source'] = 'EGM2008'
    # geoid_free2mean
    attributes['geoid_free2mean']['long_name'] = 'Free-to-Mean conversion'
    attributes['geoid_free2mean']['description'] = ('Additive value to '
        'convert geoid heights from the tide-free system to the mean-tide system')
    attributes['geoid_free2mean']['units'] = 'meters'
    attributes['geoid_free2mean']['fill_value'] = -9999.0
    attributes['geoid_free2mean']['tide_system'] = 'tide_free'

    # output variables
    dinput = {}
    dinput['lon'] = longlimit_west + np.arange(nlon)*dlon
    dinput['lat'] = latlimit_north - np.arange(nlat)*dlat
    # geoid undulation
    dinput['geoid_h'] = np.zeros((nlat,nlon), dtype=np.float32)
    geoid_h = file_contents.reshape(nlat,nlon+1)
    dinput['geoid_h'][:,:-1] = geoid_h[:,1:-1]
    # repeat values for 360
    dinput['geoid_h'][:,-1] = dinput['geoid_h'][:,0]

    # calculate Legendre polynomial of degree 2 (unnormalized)
    gridlon, gridlat = np.meshgrid(dinput['lon'],dinput['lat'])
    P2 = 0.5*(3.0*np.sin(gridlat*np.pi/180.0)**2 - 1.0)
    # offset for converting from tide_free to mean_tide
    # from Rapp 1991 (Consideration of Permanent Tidal Deformation)
    dinput['geoid_free2mean'] = np.zeros((nlat,nlon), dtype=np.float32)
    dinput['geoid_free2mean'][:,:] = -0.198*P2*(1.0 + LOVE)

    # output data and parameters to netCDF4
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    ncdf_geoid_write(dinput, attributes, FILENAME=FILENAME)
    # change permissions mode to MODE
    FILENAME.chmod(mode=MODE)

# PURPOSE: write output geoid height data to file
def ncdf_geoid_write(dinput, attributes, FILENAME=None):
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")
    # dictionary for netCDF4 variables
    nc = {}

    # defining the NetCDF dimensions
    for key in ['lon','lat']:
        fileID.createDimension(key, len(dinput[key]))
    nc['lat'] = fileID.createVariable('lat', dinput['lat'].dtype, ('lat',))
    nc['lon'] = fileID.createVariable('lon', dinput['lon'].dtype, ('lon',))
    # Defining attributes for longitude and latitude
    nc['lon'].long_name = 'longitude'
    nc['lon'].units = 'degrees_east'
    nc['lat'].long_name = 'latitude'
    nc['lat'].units = 'degrees_north'
    # defining the NetCDF functional variables
    for functional in ['geoid_h','geoid_free2mean']:
        fill_value = np.float64(attributes[functional]['fill_value'])
        nc[functional] = fileID.createVariable(functional,
            dinput[functional].dtype, ('lat','lon',),
            fill_value=fill_value, zlib=True)
        # Defining attributes for functional
        for att_name,att_val in attributes[functional].items():
            nc[functional].setncattr(att_name, att_val)

    # filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = val[:].copy()

    # Defining global attributes of NetCDF file
    for att_name,att_val in attributes['ROOT'].items():
        fileID.setncattr(key, att_val)

    # add software information
    fileID.software_reference = geoidtk.version.project_name
    fileID.software_version = geoidtk.version.full_version
    # add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()

    # Output NetCDF structure information
    logging.info(FILENAME)
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads EGM2008 geoid height spatial grids from
            unformatted binary files provided by the National
            Geospatial-Intelligence Agency
            """
    )
    # command line parameters
    parser.add_argument('gravity',
        type=pathlib.Path,
        help='Geoid height spatial grid file')
    # output filename (will default to input file with netCDF4 suffix)
    parser.add_argument('--filename','-F',
        type=pathlib.Path,
        help='Output netCDF4 filename')
    # load love number of degree 2 (default EGM2008 value)
    parser.add_argument('--love','-n',
        type=float, default=0.3,
        help='Degree 2 load Love number')
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information for each output file')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # verify input and output files
    args.gravity = pathlib.Path(args.gravity).expanduser().absolute()
    # set output file from input filename if not entered
    if not args.filename:
        args.filename = args.gravity.with_name(f'{args.gravity.stem}.nc')
    else:
        args.filename = pathlib.Path(args.filename).expanduser().absolute()

    # run program
    read_EGM2008_geoid_grids(args.gravity , FILENAME=args.filename,
        LOVE=args.love, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
