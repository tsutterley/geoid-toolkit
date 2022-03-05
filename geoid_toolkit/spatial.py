#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (03/2022)

Utilities for reading, writing and operating on spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format
        https://www.h5py.org/
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

UPDATE HISTORY:
    Updated 03/2022: add option to specify output GDAL driver
    Updated 01/2022: use iteration breaks in convert ellipsoid function
        remove fill_value attribute after creating netCDF4 and HDF5 variables
    Updated 11/2021: added empty cases to netCDF4 and HDF5 output for crs
        try to get grid mapping attributes from netCDF4 and HDF5
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: can calculate height differences between ellipsoids
    Updated 07/2021: added function for determining input variable type
    Updated 03/2021: added polar stereographic area scale calculation
        add routines for converting to and from cartesian coordinates
        eplaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: add streaming from bytes for ascii, netCDF4, HDF5, geotiff
        set default time for geotiff files to 0
    Updated 12/2020: added module for converting ellipsoids
    Updated 11/2020: output data as masked arrays if containing fill values
        add functions to read from and write to geotiff image formats
    Written 09/2020
"""
import os
import re
import io
import copy
import gzip
import uuid
import h5py
import yaml
import logging
import netCDF4
import datetime
import warnings
import numpy as np
try:
    import osgeo.gdal, osgeo.osr
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("GDAL not available")

def case_insensitive_filename(filename):
    """
    Searches a directory for a filename without case dependence
    """
    #-- check if file presently exists with input case
    if not os.access(os.path.expanduser(filename),os.F_OK):
        #-- search for filename without case dependence
        basename = os.path.basename(filename)
        directory = os.path.dirname(os.path.expanduser(filename))
        f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
        if not f:
            raise IOError('{0} not found in file system'.format(filename))
        filename = os.path.join(directory,f.pop())
    return os.path.expanduser(filename)

def data_type(x, y, t):
    """
    Determines input data type based on variable dimensions
    Inputs: spatial and temporal coordinates
    """
    xsize = np.size(x)
    ysize = np.size(y)
    tsize = np.size(t)
    if (xsize == 1) and (ysize == 1) and (tsize >= 1):
        return 'time series'
    elif (xsize == ysize) & (xsize == tsize):
        return 'drift'
    elif (np.ndim(x) > 1) & (xsize == ysize):
        return 'grid'
    elif (xsize != ysize):
        return 'grid'
    else:
        raise ValueError('Unknown data type')

def from_ascii(filename, **kwargs):
    """
    Read data from an ascii file
    Inputs: full path of input ascii file
    Options:
        ascii file is compressed or streamed from memory
        column names of ascii file
        header lines to skip from start of file
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    kwargs.setdefault('columns',['time','y','x','data'])
    kwargs.setdefault('header',0)
    #-- print filename
    logging.info(filename)
    #-- get column names
    columns = copy.copy(kwargs['columns'])
    #-- open the ascii file and extract contents
    if (kwargs['compression'] == 'gzip'):
        #-- read input ascii data from gzip compressed file and split lines
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            file_contents = f.read().decode('ISO-8859-1').splitlines()
    elif (kwargs['compression'] == 'bytes'):
        #-- read input file object and split lines
        file_contents = filename.read().splitlines()
    else:
        #-- read input ascii file (.txt, .asc) and split lines
        with open(case_insensitive_filename(filename),'r') as f:
            file_contents = f.read().splitlines()
    #-- number of lines in the file
    file_lines = len(file_contents)
    #-- compile regular expression operator for extracting numerical values
    #-- from input ascii files of spatial data
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- check if header has a known format
    if (str(kwargs['header']).upper() == 'YAML'):
        #-- counts the number of lines in the header
        YAML = False
        count = 0
        #-- Reading over header text
        while (YAML is False) & (count < file_lines):
            #-- file line at count
            line = file_contents[count]
            #-- if End of YAML Header is found: set YAML flag
            YAML = bool(re.search(r"\# End of YAML header",line))
            #-- add 1 to counter
            count += 1
        #-- parse the YAML header (specifying yaml loader)
        YAML_HEADER = yaml.load('\n'.join(file_contents[:count]),
           Loader=yaml.BaseLoader)
        #-- output spatial data and attributes
        dinput = {}
        #-- copy global attributes
        dinput['attributes'] = YAML_HEADER['header']['global_attributes']
        #-- allocate for each variable and copy variable attributes
        for c in columns:
            dinput[c] = np.zeros((file_lines-count))
            dinput['attributes'][c] = YAML_HEADER['header']['variables'][c]
        #-- update number of file lines to skip for reading data
        header = int(count)
    else:
        #-- output spatial data and attributes
        dinput = {c:np.zeros((file_lines-kwargs['header'])) for c in columns}
        dinput['attributes'] = {c:dict() for c in columns}
    #-- extract spatial data array
    #-- for each line in the file
    for i,line in enumerate(file_contents[header:]):
        #-- extract columns of interest and assign to dict
        #-- convert fortran exponentials if applicable
        column = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
        #-- copy variables from column dict to output dictionary
        for c in columns:
            dinput[c][i] = np.float64(column[c])
    #-- convert to masked array if fill values
    dinput['data'] = np.ma.asarray(dinput['data'])
    dinput['data'].mask = np.zeros_like(dinput['data'],dtype=bool)
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask[:] = (dinput['data'].data == dinput['data'].fill_value)
    #-- return the spatial variables
    return dinput

def from_netCDF4(filename, **kwargs):
    """
    Read data from a netCDF4 file
    Inputs: full path of input netCDF4 file
    Options:
        netCDF4 file is compressed or streamed from memory
        netCDF4 variable names of time, longitude, latitude, and data
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    kwargs.setdefault('timename','time')
    kwargs.setdefault('xname','lon')
    kwargs.setdefault('yname','lat')
    kwargs.setdefault('varname','data')
    #-- read data from netCDF4 file
    #-- Open the NetCDF4 file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=f.read())
    elif (kwargs['compression'] == 'bytes'):
        #-- read as in-memory (diskless) netCDF4 dataset
        fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=filename.read())
    else:
        #-- read netCDF4 dataset
        fileID = netCDF4.Dataset(case_insensitive_filename(filename), 'r')
    #-- Output NetCDF file information
    logging.info(fileID.filepath())
    logging.info(list(fileID.variables.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            ncattr, = [s for s in fileID.ncattrs() if re.match(attr,s,re.I)]
            dinput['attributes'][attr] = fileID.getncattr(ncattr)
        except (ValueError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','grid_mapping','_FillValue']
    #-- mapping between netCDF4 variable names and output names
    variable_mapping = dict(x=kwargs['xname'],y=kwargs['yname'],
        data=kwargs['varname'],time=kwargs['timename'])
    #-- for each variable
    for key,nc in variable_mapping.items():
        #-- Getting the data from each NetCDF variable
        dinput[key] = fileID.variables[nc][:]
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                ncattr, = [s for s in fileID.variables[nc].ncattrs()
                    if re.match(attr,s,re.I)]
                dinput['attributes'][key][attr] = \
                    fileID.variables[nc].getncattr(ncattr)
            except (ValueError,AttributeError):
                pass
    #-- get projection information if there is a grid_mapping attribute
    if 'grid_mapping' in dinput['attributes']['data'].keys():
        #-- try getting the attribute
        grid_mapping = dinput['attributes']['data']['grid_mapping']
        for att_name in fileID[grid_mapping].ncattrs():
            dinput['attributes']['crs'][att_name] = \
                fileID.variables[grid_mapping].getncattr(att_name)
        #-- get the spatial projection reference information from wkt
        #-- and overwrite the file-level projection attribute (if existing)
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(dinput['attributes']['crs']['crs_wkt'])
        dinput['attributes']['projection'] = srs.ExportToProj4()
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- Closing the NetCDF file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_HDF5(filename, **kwargs):
    """
    Read data from a HDF5 file
    Inputs: full path of input HDF5 file
    Options:
        HDF5 file is compressed or streamed from memory
        HDF5 variable names of time, longitude, latitude, and data
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    kwargs.setdefault('timename','time')
    kwargs.setdefault('xname','lon')
    kwargs.setdefault('yname','lat')
    kwargs.setdefault('varname','data')
    #-- read data from HDF5 file
    #-- Open the HDF5 file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read gzip compressed file and extract into in-memory file object
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fid = io.BytesIO(f.read())
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (kwargs['compression'] == 'bytes'):
        #-- read as in-memory (diskless) HDF5 dataset
        fileID = h5py.File(filename, 'r')
    else:
        #-- read HDF5 dataset
        fileID = h5py.File(case_insensitive_filename(filename), 'r')
    #-- Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            dinput['attributes'][attr] = fileID.attrs[attr]
        except (KeyError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','grid_mapping','_FillValue']
    #-- mapping between HDF5 variable names and output names
    variable_mapping = dict(x=kwargs['xname'],y=kwargs['yname'],
        data=kwargs['varname'],time=kwargs['timename'])
    #-- for each variable
    for key,h5 in variable_mapping.items():
        #-- Getting the data from each HDF5 variable
        dinput[key] = np.copy(fileID[h5][:])
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                dinput['attributes'][key][attr] = fileID[h5].attrs[attr]
            except (KeyError,AttributeError):
                pass
    #-- get projection information if there is a grid_mapping attribute
    if 'grid_mapping' in dinput['attributes']['data'].keys():
        #-- try getting the attribute
        grid_mapping = dinput['attributes']['data']['grid_mapping']
        for att_name,att_val in fileID[grid_mapping].attrs.items():
            dinput['attributes']['crs'][att_name] = att_val
        #-- get the spatial projection reference information from wkt
        #-- and overwrite the file-level projection attribute (if existing)
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(dinput['attributes']['crs']['crs_wkt'])
        dinput['attributes']['projection'] = srs.ExportToProj4()
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- Closing the HDF5 file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_geotiff(filename, **kwargs):
    """
    Read data from a geotiff file
    Inputs: full path of input geotiff file
    Options:
        geotiff file is compressed or streamed from memory
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    #-- Open the geotiff file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read gzip compressed file and extract into memory-mapped object
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            osgeo.gdal.FileFromMemBuffer(mmap_name, f.read())
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        ds = osgeo.gdal.Open(mmap_name)
    elif (kwargs['compression'] == 'bytes'):
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        osgeo.gdal.FileFromMemBuffer(mmap_name, filename.read())
        ds = osgeo.gdal.Open(mmap_name)
    else:
        #-- read geotiff dataset
        ds = osgeo.gdal.Open(case_insensitive_filename(filename))
    #-- print geotiff file if verbose
    logging.info(filename)
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {c:dict() for c in ['x','y','data']}
    #-- get the spatial projection reference information
    srs = ds.GetSpatialRef()
    dinput['attributes']['projection'] = srs.ExportToProj4()
    dinput['attributes']['wkt'] = srs.ExportToWkt()
    #-- get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    bsize = ds.RasterCount
    #-- get geotiff info
    info_geotiff = ds.GetGeoTransform()
    dinput['attributes']['spacing'] = (info_geotiff[1],info_geotiff[5])
    #-- calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    dinput['attributes']['extent'] = (xmin,xmax,ymin,ymax)
    #-- x and y pixel center coordinates (converted from upper left)
    dinput['x'] = xmin + info_geotiff[1]/2.0 + np.arange(xsize)*info_geotiff[1]
    dinput['y'] = ymax + info_geotiff[5]/2.0 + np.arange(ysize)*info_geotiff[5]
    #-- read full image with GDAL
    dinput['data'] = ds.ReadAsArray()
    #-- set default time to zero for each band
    dinput.setdefault('time', np.zeros((bsize)))
    #-- check if image has fill values
    dinput['data'] = np.ma.asarray(dinput['data'])
    dinput['data'].mask = np.zeros_like(dinput['data'],dtype=bool)
    if ds.GetRasterBand(1).GetNoDataValue():
        #-- mask invalid values
        dinput['data'].fill_value = ds.GetRasterBand(1).GetNoDataValue()
        #-- create mask array for bad values
        dinput['data'].mask[:] = (dinput['data'].data == dinput['data'].fill_value)
        #-- set attribute for fill value
        dinput['attributes']['data']['_FillValue'] = dinput['data'].fill_value
    #-- close the dataset
    ds = None
    #-- return the spatial variables
    return dinput

def to_ascii(output, attributes, filename, **kwargs):
    """
    Write data to an ascii file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output ascii file
    Options:
        delimiter for output spatial file
        order of columns for output spatial file
        create a YAML header with data attributes
    """
    #-- set default keyword arguments
    kwargs.setdefault('delimiter',',')
    kwargs.setdefault('columns',['time','lat','lon','tide'])
    kwargs.setdefault('header',False)
    #-- get column names
    columns = copy.copy(kwargs['columns'])
    #-- output filename
    filename = os.path.expanduser(filename)
    logging.info(filename)
    #-- open the output file
    fid = open(filename, 'w')
    #-- create a column stack arranging data in column order
    data_stack = np.c_[[output[col] for col in columns]]
    ncol,nrow = np.shape(data_stack)
    #-- print YAML header to top of file
    if kwargs['header']:
        fid.write('{0}:\n'.format('header'))
        #-- data dimensions
        fid.write('\n  {0}:\n'.format('dimensions'))
        fid.write('    {0:22}: {1:d}\n'.format('time',nrow))
        #-- non-standard attributes
        fid.write('  {0}:\n'.format('non-standard_attributes'))
        #-- data format
        fid.write('    {0:22}: ({1:d}f0.8)\n'.format('formatting_string',ncol))
        fid.write('\n')
        #-- global attributes
        fid.write('\n  {0}:\n'.format('global_attributes'))
        today = datetime.datetime.now().isoformat()
        fid.write('    {0:22}: {1}\n'.format('date_created', today))
        # print variable descriptions to YAML header
        fid.write('\n  {0}:\n'.format('variables'))
        #-- print YAML header with variable attributes
        for i,v in enumerate(columns):
            fid.write('    {0:22}:\n'.format(v))
            for atn,atv in attributes[v].items():
                fid.write('      {0:20}: {1}\n'.format(atn,atv))
            #-- add precision and column attributes for ascii yaml header
            fid.write('      {0:20}: double_precision\n'.format('precision'))
            fid.write('      {0:20}: column {1:d}\n'.format('comment',i+1))
        #-- end of header
        fid.write('\n\n# End of YAML header\n')
    #-- write to file for each data point
    for line in range(nrow):
        line_contents = ['{0:0.8f}'.format(d) for d in data_stack[:,line]]
        print(kwargs['delimiter'].join(line_contents), file=fid)
    #-- close the output file
    fid.close()

def to_netCDF4(output, attributes, filename, **kwargs):
    """
    Write data to a netCDF4 file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output netCDF4 file
    """
    #-- opening NetCDF file for writing
    fileID = netCDF4.Dataset(os.path.expanduser(filename),'w',format="NETCDF4")
    #-- Defining the NetCDF dimensions
    fileID.createDimension('time', len(np.atleast_1d(output['time'])))
    #-- defining the NetCDF variables
    nc = {}
    for key,val in output.items():
        if '_FillValue' in attributes[key].keys():
            nc[key] = fileID.createVariable(key, val.dtype, ('time',),
                fill_value=attributes[key]['_FillValue'], zlib=True)
            attributes[key].pop('_FillValue')
        elif val.shape:
            nc[key] = fileID.createVariable(key, val.dtype, ('time',))
        else:
            nc[key] = fileID.createVariable(key, val.dtype, ())
        #-- filling NetCDF variables
        nc[key][:] = val
        #-- Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            nc[key].setncattr(att_name,att_val)
    #-- add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()
    #-- add file-level attributes if applicable
    if 'ROOT' in attributes.keys():
        #-- Defining attributes for file
        for att_name,att_val in attributes['ROOT'].items():
            fileID.setncattr(att_name,att_val)
    #-- Output NetCDF structure information
    logging.info(filename)
    logging.info(list(fileID.variables.keys()))
    #-- Closing the NetCDF file
    fileID.close()

def to_HDF5(output, attributes, filename, **kwargs):
    """
    Write data to a HDF5 file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output HDF5 file
    """
    #-- opening HDF5 file for writing
    fileID = h5py.File(filename, 'w')
    #-- Defining the HDF5 dataset variables
    h5 = {}
    for key,val in output.items():
        if '_FillValue' in attributes[key].keys():
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, fillvalue=attributes[key]['_FillValue'],
                compression='gzip')
            attributes[key].pop('_FillValue')
        elif val.shape:
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, compression='gzip')
        else:
            h5[key] = fileID.create_dataset(key, val.shape,
                dtype=val.dtype)
        #-- Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            h5[key].attrs[att_name] = att_val
    #-- add attribute for date created
    fileID.attrs['date_created'] = datetime.datetime.now().isoformat()
    #-- add file-level attributes if applicable
    if 'ROOT' in attributes.keys():
        #-- Defining attributes for file
        for att_name,att_val in attributes['ROOT'].items():
            fileID.attrs[att_name] = att_val
    #-- Output HDF5 structure information
    logging.info(filename)
    logging.info(list(fileID.keys()))
    #-- Closing the HDF5 file
    fileID.close()

def to_geotiff(output, attributes, filename, **kwargs):
    """
    Write data to a geotiff file
    Inputs:
        python dictionary of output data
        python dictionary of output attributes
        full path of output geotiff file
    Options:
        output variable name
        GDAL driver
        GDAL data type
        GDAL driver creation options
    """
    #-- set default keyword arguments
    kwargs.setdefault('varname','data')
    kwargs.setdefault('driver',"GTiff")
    kwargs.setdefault('dtype',osgeo.gdal.GDT_Float64)
    kwargs.setdefault('options',['COMPRESS=LZW'])
    varname = copy.copy(kwargs['varname'])
    #-- verify grid dimensions to be iterable
    output = expand_dims(output, varname=varname)
    #-- grid shape
    ny,nx,nband = np.shape(output[varname])
    #-- output as geotiff or specified driver
    driver = osgeo.gdal.GetDriverByName(kwargs['driver'])
    #-- set up the dataset with creation options
    ds = driver.Create(filename, nx, ny, nband,
        kwargs['dtype'], kwargs['options'])
    #-- top left x, w-e pixel resolution, rotation
    #-- top left y, rotation, n-s pixel resolution
    xmin,xmax,ymin,ymax = attributes['extent']
    dx,dy = attributes['spacing']
    ds.SetGeoTransform([xmin,dx,0,ymax,0,dy])
    #-- set the spatial projection reference information
    srs = osgeo.osr.SpatialReference()
    srs.ImportFromWkt(attributes['wkt'])
    #-- export
    ds.SetProjection( srs.ExportToWkt() )
    #-- for each band
    for band in range(nband):
        #-- set fill value for band
        if '_FillValue' in attributes[varname].keys():
            fill_value = attributes[varname]['_FillValue']
            ds.GetRasterBand(band+1).SetNoDataValue(fill_value)
        #-- write band to geotiff array
        ds.GetRasterBand(band+1).WriteArray(output[varname][:,:,band])
    #-- print filename if verbose
    logging.info(filename)
    #-- close dataset
    ds.FlushCache()

def expand_dims(obj, varname='data'):
    """
    Add a singleton dimension to a spatial dictionary if non-existent
    Options:
        variable name to modify
    """
    #-- change time dimensions to be iterableinformation
    try:
        obj['time'] = np.atleast_1d(obj['time'])
    except:
        pass
    #-- output spatial with a third dimension
    if isinstance(varname,list):
        for v in varname:
            obj[v] = np.atleast_3d(obj[v])
    elif isinstance(varname,str):
        obj[varname] = np.atleast_3d(obj[varname])
    #-- return reformed spatial dictionary
    return obj

def convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10):
    """
    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Inputs:
        phi1: latitude of input ellipsoid in degrees
        h1: height above input ellipsoid in meters
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid

    Options:
        eps: tolerance to prevent division by small numbers
            and to determine convergence
        itmax: maximum number of iterations to use in Newton-Raphson

    Returns:
        phi2: latitude of output ellipsoid in degrees
        h2: height above output ellipsoid in meters

    References:
        Astronomical Algorithms, Jean Meeus, 1991, Willmann-Bell, Inc.
            pp. 77-82
    """
    if (len(phi1) != len(h1)):
        raise ValueError('phi and h have incompatable dimensions')
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- initialize output arrays
    npts = len(phi1)
    phi2 = np.zeros((npts))
    h2 = np.zeros((npts))
    #-- for each point
    for N in range(npts):
        #-- force phi1 into range -90 <= phi1 <= 90
        if (np.abs(phi1[N]) > 90.0):
            phi1[N] = np.sign(phi1[N])*90.0
        #-- handle special case near the equator
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + a1 - a2
        if (np.abs(phi1[N]) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + a1 - a2
        #-- handle special case near the poles
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + b1 - b2
        elif ((90.0 - np.abs(phi1[N])) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + b1 - b2
        #-- handle case if latitude is within 45 degrees of equator
        elif (np.abs(phi1[N]) <= 45):
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = b2 * b2 - a2 * a2
            k1 = a2 * hpr1cos
            k2 = b2 * hpr1sin
            #-- perform newton-raphson iteration to solve for u2
            #-- cos(u2) will not be close to zero since abs(phi1) <= 45
            for i in range(0, itmax+1):
                cosu2 = np.cos(u2)
                fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                if (np.abs(fu2p) < eps):
                    break
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        break
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(phi2r)
        #-- handle final case where latitudes are between 45 degrees and pole
        else:
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = a2 * a2 - b2 * b2
            k1 = b2 * hpr1sin
            k2 = a2 * hpr1cos
            #-- perform newton-raphson iteration to solve for u2
            #-- sin(u2) will not be close to zero since abs(phi1) > 45
            for i in range(0, itmax+1):
                sinu2 = np.sin(u2)
                fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                fu2p =  -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                if (np.abs(fu2p) < eps):
                    break
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        break
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(phi2r)

    #-- return the latitude and height
    return (phi2, h2)

def compute_delta_h(a1, f1, a2, f2, lat):
    """
    Compute difference in elevation for two ellipsoids at a given
        latitude using a simplified empirical equation

    Inputs:
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid
        lat: array of latitudes in degrees

    Returns:
        delta_h: difference in elevation for two ellipsoids

    Reference:
        J Meeus, Astronomical Algorithms, pp. 77-82 (1991)
    """
    #-- force phi into range -90 <= phi <= 90
    gt90, = np.nonzero((lat < -90.0) | (lat > 90.0))
    lat[gt90] = np.sign(lat[gt90])*90.0
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- compute delta_a and delta_b coefficients
    delta_a = a2 - a1
    delta_b = b2 - b1
    #-- compute differences between ellipsoids
    #-- delta_h = -(delta_a * cos(phi)^2 + delta_b * sin(phi)^2)
    phi = lat * np.pi/180.0
    delta_h = -(delta_a*np.cos(phi)**2 + delta_b*np.sin(phi)**2)
    return delta_h

def wrap_longitudes(lon):
    """
    Wraps longitudes to range from -180 to +180

    Inputs:
        lon: longitude (degrees east)
    """
    phi = np.arctan2(np.sin(lon*np.pi/180.0),np.cos(lon*np.pi/180.0))
    #-- convert phi from radians to degrees
    return phi*180.0/np.pi

def to_cartesian(lon,lat,h=0.0,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Converts geodetic coordinates to Cartesian coordinates

    Inputs:
        lon: longitude (degrees east)
        lat: latitude (degrees north)

    Options:
        h: height above ellipsoid (or sphere)
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
            * for spherical coordinates set to radius of the Earth
        flat: ellipsoidal flattening (default: WGS84)
            * for spherical coordinates set to 0
    """
    #-- verify axes
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    #-- fix coordinates to be 0:360
    lon[lon < 0] += 360.0
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- convert from geodetic latitude to geocentric latitude
    dtr = np.pi/180.0
    #-- geodetic latitude in radians
    latitude_geodetic_rad = lat*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(lon*dtr)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(lon*dtr)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    #-- return the cartesian coordinates
    return (X,Y,Z)

def to_sphere(x,y,z):
    """
    Convert from cartesian coordinates to spherical coordinates

    Inputs:
        x,y,z in cartesian coordinates
    """
    #-- calculate radius
    rad = np.sqrt(x**2.0 + y**2.0 + z**2.0)
    #-- calculate angular coordinates
    #-- phi: azimuthal angle
    phi = np.arctan2(y,x)
    #-- th: polar angle
    th = np.arccos(z/rad)
    #-- convert to degrees and fix to 0:360
    lon = 180.0*phi/np.pi
    if np.any(lon < 0):
        lt0 = np.nonzero(lon < 0)
        lon[lt0] += 360.0
    #-- convert to degrees and fix to -90:90
    lat = 90.0 - (180.0*th/np.pi)
    np.clip(lat, -90, 90, out=lat)
    #-- return latitude, longitude and radius
    return (lon,lat,rad)

def to_geodetic(x,y,z,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using a closed form solution

    Inputs:
        x,y,z in cartesian coordinates

    Options:
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
        flat: ellipsoidal flattening (default: WGS84)

    References:
        J Zhu "Exact conversion of Earth-centered, Earth-fixed
            coordinates to geodetic coordinates"
        Journal of Guidance, Control, and Dynamics,
        16(2), 389--391, 1993
        https://arc.aiaa.org/doi/abs/10.2514/3.21016
    """
    #-- semiminor axis of the WGS84 ellipsoid [m]
    b_axis = (1.0 - flat)*a_axis
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- square of first numerical eccentricity
    e12 = ecc1**2
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- calculate distance
    w = np.sqrt(x**2 + y**2)
    #-- calculate longitude
    lon = np.arctan2(y,x)/dtr
    lat = np.zeros_like(lon)
    h = np.zeros_like(lon)
    if (w == 0):
        #-- special case where w == 0 (exact polar solution)
        h = np.sign(z)*z - b_axis
        lat = 90.0*np.sign(z)
    else:
        #-- all other cases
        l = e12/2.0
        m = (w/a_axis)**2.0
        n = ((1.0-e12)*z/b_axis)**2.0
        i = -(2.0*l**2 + m + n)/2.0
        k = (l**2.0 - m - n)*l**2.0
        q = (1.0/216.0)*(m + n - 4.0*l**2)**3.0 + m*n*l**2.0
        D = np.sqrt((2.0*q - m*n*l**2)*m*n*l**2)
        B = i/3.0 - (q+D)**(1.0/3.0) - (q-D)**(1.0/3.0)
        t = np.sqrt(np.sqrt(B**2-k) - (B+i)/2.0)-np.sign(m-n)*np.sqrt((B-i)/2.0)
        wi = w/(t+l)
        zi = (1.0-e12)*z/(t-l)
        #-- calculate latitude and height
        lat = np.arctan2(zi,((1.0-e12)*wi))/dtr
        h = np.sign(t-1.0+l)*np.sqrt((w-wi)**2.0 + (z-zi)**2.0)
    #-- return latitude, longitude and height
    return (lon,lat,h)
