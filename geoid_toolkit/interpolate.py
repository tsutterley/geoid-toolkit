"""
interpolate.py
Written by Tyler Sutterley (07/2026)
Routines to interpolate data from pre-computed spatial grids

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/

UPDATE HISTORY:
    Written 07/2026
"""

from __future__ import annotations

import logging
import numpy as np
from geoid_toolkit.utilities import import_dependency, get_data_path

# attempt imports
netCDF4 = import_dependency('netCDF4')

__all__ = [
    'Interpolate',
    'geoid_height',
]


class Interpolate:
    """
    Class for interpolating data from pre-computed spatial grids
    """

    def __init__(self, coords, values=None, **kwargs):
        self.xp = np.copy(coords[1])
        self.yp = np.copy(coords[0])
        self.shape = (len(self.yp), len(self.xp))
        self.fp = np.zeros(self.shape)
        # default interpolation and extrapolation method
        self.method = 'linear'
        self.extrapolate = 'linear'
        self.fill_value = np.nan
        # update values if provided
        if values is not None:
            self.update(values)
        # set any keyword arguments
        for key, val in kwargs.items():
            setattr(self, key, val)

    def update(self, values):
        """Update the interpolation grid values"""
        if getattr(self, '_flip_y', False):
            self.fp[:] = values[::-1, :]
        else:
            self.fp[:] = np.copy(values)
        # validate the shape of the input values
        self.__validate__()

    def __call__(self, points, **kwargs):
        """Run interpolation for input points and values

        Parameters
        ----------
        points: tuple
            (y, x) tuple of arrays containing the points to interpolate
        """
        if np.size(points[1]) != np.size(points[0]):
            # assume points are 2D arrays and create meshgrid
            y, x = np.meshgrid(points[0], points[1], indexing='ij')
        else:
            # extract points
            y, x = np.copy(points)
        # set default keyword arguments
        fp = kwargs.get('values', self.fp)
        kwargs.setdefault('method', self.method)
        kwargs.setdefault('extrapolate', self.extrapolate)
        kwargs.setdefault('fill_value', self.fill_value)
        # run interpolation
        return self._interp(
            (y.flatten(), x.flatten()),
            (self.yp, self.xp),
            fp,
            **kwargs,
        ).reshape(y.shape)

    def __validate__(self):
        """Validate the shape of the input values"""
        if self.fp.shape != self.shape:
            raise ValueError('Input shape does not match expected')
        # set longitudinal convention to be 0:360 degrees
        if self.xp.min() < 0.0:
            self.xp = np.where(self.xp < 0.0, self.xp + 360.0, self.xp)
        # check that the grid is increasing in the y-direction
        if self.yp[1] - self.yp[0] <= 0.0:
            self.yp = self.yp[::-1]
            self.fp = self.fp[::-1, :]
            # set flag to flip over y-coordinates for future updates
            self._flip_y = True

    @staticmethod
    def _interp(points, coords, fp, **kwargs):
        """
        Interpolate data from pre-computed spatial grids

        Parameters
        ----------
        points: tuple
            Tuple of arrays containing the points to interpolate
        coords: tuple
            Tuple of arrays containing the coordinates of the grid
        fp: np.ndarray
            Array of values on the grid to interpolate
        method: str, default "linear"
            Interpolation method to use

                - ``'linear'``: linear interpolation
                - ``'nearest'``: nearest-neighbor interpolation
        extrapolate: str, default "linear"
            Extrapolation method to use

                - ``'linear'``: linearly extrapolate
                - ``'nearest'``: use nearest-neighbor value
                - ``'fill_value'``: replace with fill value
                - ``'zero'``: replace with zero
        fill_value: float, default np.nan
            Invalid values

        Returns
        -------
        f: np.ndarray
            Interpolated values
        """
        # extract points and coordinates
        y, x = points
        yp, xp = coords
        # get interpolation method
        method = kwargs.get('method', 'linear').lower()
        if method not in ('linear', 'nearest'):
            raise ValueError(f'Invalid interpolation method: {method}')
        # get extrapolation method
        extrapolate = kwargs.get('extrapolate', 'linear').lower()
        if extrapolate not in ('linear', 'nearest', 'fill_value', 'zero'):
            raise ValueError(f'Invalid extrapolate method: {extrapolate}')
        # replace fill values with NaN for interpolation
        fill_value = kwargs.get('fill_value', np.nan)
        fp = np.where((fp == fill_value) | np.isnan(fp), np.nan, fp)
        # grid spacing in longitudinal direction
        dx = np.abs(xp[1] - xp[0])
        # adjust x-coordinates to be consistent with model
        if (np.min(x) < 0.0) & (xp.max() > (180.0 + dx)):
            # input points convention (-180:180)
            # geoid model convention (0:360)
            x = np.where(x < 0.0, x + 360.0, x)
        # find indices outside of the grid for extrapolation
        if extrapolate in ('zero', 'fill_value'):
            # find indices outside of the grid
            outside = (
                (y < yp.min())
                | (y > yp.max())
                | (x < xp.min())
                | (x > xp.max())
            )
        # clip coordinates to handle nearest-neighbor extrapolation
        if extrapolate in ('nearest', 'zero', 'fill_value'):
            y = np.clip(y, a_min=yp.min(), a_max=yp.max())
            x = np.clip(x, a_min=xp.min(), a_max=xp.max())
        # find indice where y/x could be inserted into yp/xp
        i = np.searchsorted(yp, y) - 1
        j = np.searchsorted(xp, x) - 1
        # clip indices to handle linear extrapolation
        if extrapolate == 'linear':
            i = np.clip(i, a_min=0, a_max=len(yp) - 2)
            j = np.clip(j, a_min=0, a_max=len(xp) - 2)
        # fractional distance between points
        disty = np.divide(y - yp[i], yp[i + 1] - yp[i])
        distx = np.divide(x - xp[j], xp[j + 1] - xp[j])
        # round if using nearest-neighbors method
        if method == 'nearest':
            disty = np.around(disty).astype(np.int32)
            distx = np.around(distx).astype(np.int32)
        # calculate interpolated values
        f = (
            ((1.0 - disty) * (1.0 - distx) * fp[i, j])
            + ((1.0 - disty) * distx * fp[i, j + 1])
            + (disty * (1.0 - distx) * fp[i + 1, j])
            + (disty * distx * fp[i + 1, j + 1])
        )
        # replace NaN values with fill_value
        if (fill_value is not None) and not np.isnan(fill_value):
            f = np.where(np.isnan(f), fill_value, f)
        # replace values outside of the grid with fill_value or zero
        if extrapolate == 'zero':
            f = np.where(outside, 0.0, f)
        elif extrapolate == 'fill_value':
            f = np.where(outside, fill_value, f)
        # return the interpolated values
        return f


# PURPOSE: read geoid height netCDF4 file and interpolate to points
def geoid_height(
    longitude: np.ndarray,
    latitude: np.ndarray,
    model: str = 'EGM2008',
    tide_system: str = 'mean_tide',
    method: str = 'linear',
    **kwargs,
):
    """
    Interpolate geoid heights to input lat/lon points

    Parameters
    ----------
    longitude: np.ndarray
        Longitudes of output points
    latitude: np.ndarray
        Latitudes of output points
    model: str, default "EGM2008"
        Geoid model to use in the interpolation
    tide_system: str, default "mean_tide"
        Output tide-system for the geoid model

        - ``'mean_tide'``: mean-tide system
        - ``'tide_free'``: tide-free system
        - ``'zero_tide'``: zero-tide system
    method: str, default "linear"
        Interpolation method to use
        
        - ``'linear'``: linear interpolation
        - ``'nearest'``: nearest-neighbor interpolation
    kwargs: keyword arguments
        Additional keyword arguments for the ``Interpolate`` class

    Returns
    -------
    geoid_undulation: np.ndarray
        Geoid undulation at input lat/lon points
    """
    # read geoid heights and free-to-mean conversion
    FILENAME = get_data_path(['data', f'{model}_geoid_h.nc'])
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        geoid_h = fileID.variables['geoid_h'][:].copy()
        geoid_free2mean = fileID.variables['geoid_free2mean'][:].copy()
        lon = fileID.variables['lon'][:].copy()
        lat = fileID.variables['lat'][:].copy()
        # get degree 2 load Love number for free-to-mean conversion
        # assume a default value of 0.3 if not present in the file
        k2 = getattr(fileID.variables['geoid_free2mean'], 'k2', 0.3)
    # interpolate values for a given method
    interp = Interpolate((lat, lon), geoid_h, method=method, **kwargs)
    # interpolate geoid heights for input lat/lon points
    geoid_undulation = interp((latitude, longitude))
    if tide_system.lower() == 'mean_tide':
        # update values for interpolator
        interp.update(geoid_free2mean)
        # convert to mean-tide system using free-to-mean conversion
        geoid_undulation += interp((latitude, longitude))
    elif tide_system.lower() == 'zero_tide':
        # update values for interpolator
        interp.update(geoid_free2mean)
        # convert to zero-tide system using free-to-mean conversion and k2
        geoid_undulation += interp((latitude, longitude)) / (1.0 + k2)
    # return geoid undulation
    return geoid_undulation
