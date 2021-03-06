#!/usr/bin/env python
u"""
read_ICGEM_geoid_grids.py
Written by Tyler Sutterley (03/2021)
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
import os
import re
import netCDF4
import argparse
import numpy as np

#-- PURPOSE: Reads .gdf grids from the GFZ calculation service
def read_ICGEM_geoid_grids(FILE, FILENAME=None, MARKER='', SPACING=None,
    VERBOSE=False, MODE=0o775):
    #-- split filename into basename and extension
    fileBasename,_ = os.path.splitext(FILE)
    #-- open input file and read contents
    with open(FILE, 'r') as f:
        file_contents = f.read().splitlines()
    #-- number of lines contained in the file
    file_lines = len(file_contents)

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text and extracting parameters
    HEADER = False
    parameters = {}
    while HEADER is False:
        #-- file line at count
        line = file_contents[count]
        if (len(line) > 1):
            col = line.split()
            parameters[col[0].strip()] = col[1].strip()
        #-- find MARKER within line to set HEADER flag to True when found
        HEADER = bool(re.match(MARKER,line))
        #-- add 1 to counter
        count += 1

    #-- clean up dictionary of parameters
    for key in ['[deg.]','longitude',MARKER]:
        parameters = removekey(parameters, key)

    #-- extract necessary parameters
    latlimit_north = np.float(parameters['latlimit_north'])
    latlimit_south = np.float(parameters['latlimit_south'])
    longlimit_west = np.float(parameters['longlimit_west'])
    longlimit_east = np.float(parameters['longlimit_east'])
    #-- change grid spacing by binning data
    if SPACING is None:
        nlat = np.int(parameters['latitude_parallels'])
        nlon = np.int(parameters['longitude_parallels'])
        dlon = np.float(parameters['gridstep'])
        dlat = np.float(parameters['gridstep'])
    else:
        dlon,dlat = SPACING
        parameters['gridstep'] = '{0:g},{1:g}'.format(dlon,dlat)
        nlat = np.abs((latlimit_north-latlimit_south)/dlat).astype('i') + 1
        nlon = np.abs((longlimit_west-longlimit_east)/dlon).astype('i') + 1
        parameters['latitude_parallels'] = '{0:d}'.format(nlat)
        parameters['longitude_parallels'] = '{0:d}'.format(nlon)

    #-- output dataset
    dinput = {}
    functional = parameters['functional']
    dinput[functional] = np.zeros((nlat,nlon))
    dinput['lon'] = longlimit_west + np.arange(nlon)*dlon
    dinput['lat'] = latlimit_north - np.arange(nlat)*dlat

    #-- for each file line
    bin_count = np.zeros((nlat,nlon))
    for j in range(count, file_lines):
        col = np.array(file_contents[j].split(), dtype=np.float)
        #-- calculating the lon/lat indice
        ilon = np.int((longlimit_west + col[0])/dlon)
        ilat = np.int((latlimit_north - col[1])/dlat)
        #-- if wanting data lat/lon
        dinput[functional][ilat,ilon] += np.float(col[2])
        bin_count[ilat,ilon] += 1.0

    #-- take the mean of the binned data (if not regridding will divide by 1)
    ii,jj = np.nonzero(bin_count > 0)
    dinput[functional][ii,jj] /= bin_count[ii,jj]
    ii,jj = np.nonzero(bin_count == 0)
    dinput[functional][ii,jj] = np.float(parameters['gapvalue'])

    #-- output data and parameters to netCDF4
    FILENAME = '{0}.nc'.format(fileBasename) if (FILENAME is None) else FILENAME
    ncdf_geoid_write(dinput, parameters, FILENAME=FILENAME, VERBOSE=VERBOSE)
    #-- change permissions mode to MODE
    os.chmod(FILENAME, MODE)

#-- PURPOSE: remove keys from an input dictionary
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

#-- PURPOSE: write output geoid height data to file
def ncdf_geoid_write(dinput, parameters, FILENAME=None, VERBOSE=False):
    #-- opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    #-- Defining the NetCDF dimensions
    for key in ['lon','lat']:
        fileID.createDimension(key, len(dinput[key]))

    #-- defining the NetCDF variables
    nc = {}
    functional = parameters['functional']
    gapvalue = np.float(parameters['gapvalue'])
    nc['lat'] = fileID.createVariable('lat', dinput['lat'].dtype, ('lat',))
    nc['lon'] = fileID.createVariable('lon', dinput['lon'].dtype, ('lon',))
    nc[functional] = fileID.createVariable(functional, dinput[functional].dtype,
        ('lat','lon',), fill_value=gapvalue, zlib=True)
    #-- filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = val[:].copy()

    #-- Defining attributes for longitude and latitude
    nc['lon'].long_name = 'longitude'
    nc['lon'].units = 'degrees_east'
    nc['lat'].long_name = 'latitude'
    nc['lat'].units = 'degrees_north'
    #-- Defining attributes for functional
    nc[functional].units = parameters['unit']
    #-- global variables of NetCDF file
    for key in sorted(parameters.keys()):
        fileID.setncattr(key, parameters[key])

    #-- Output NetCDF structure information
    if VERBOSE:
        print(os.path.basename(FILENAME))
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()

#-- Main program that calls read_ICGEM_geoid_grids()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads geoid height spatial grids from the ICGEM
            Geoid Calculation Service
            """
    )
    #-- command line parameters
    parser.add_argument('files',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Geoid height spatial grid files (*.gdf)')
    #-- output filename (will default to input file with netCDF4 suffix)
    parser.add_argument('--filename','-F',
        type=str, default=None,
        help='Output netCDF4 filename')
    #-- marker denoting the end of the header text
    parser.add_argument('--header','-H',
        type=str, default='end_of_head',
        help='Marker denoting the end of the header text')
    #-- change the output grid spacing by binning
    parser.add_argument('--spacing','-S',
        type=float, default=None, nargs=2,
        help='Change output grid spacing')
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information for each output file')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    args = parser.parse_args()

    #-- for each input grid file
    for f in args.files:
        read_ICGEM_geoid_grids(f, FILENAME=args.filename, MARKER=args.header,
            SPACING=args.spacing, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
