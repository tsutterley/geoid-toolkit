"""
test_interpolate.py
"""
import pytest
import netCDF4
import numpy as np
import scipy.interpolate
import geoid_toolkit as geoidtk


# parametrize over interpolation methods
@pytest.mark.parametrize('method', ['linear', 'nearest'])
def test_interpolate_geoid(method, model='EGM2008'):
    # create a simple grid
    dlon, dlat = 1.2, 1.2
    ilat = 90.0 - dlat * (np.arange(180 // dlat) + 0.5)
    ilon = dlon * (np.arange(360 // dlon) + 0.5)
    gridlon, gridlat = np.meshgrid(ilon, ilat)
    # interpolate geoid height to grid
    test = geoidtk.interpolate.geoid_height(
        ilon, ilat, model=model, method=method, tide_system='tide_free'
    )
    # read geoid heights and free-to-mean conversion
    FILENAME = geoidtk.utilities.get_data_path(['data', f'{model}_geoid_h.nc'])
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        geoid_h = fileID.variables['geoid_h'][:].copy()
        lon = fileID.variables['lon'][:].copy()
        lat = fileID.variables['lat'][:].copy()
    # interpolate using regular grid interpolation
    rgi = scipy.interpolate.RegularGridInterpolator(
        (lat, lon), geoid_h, method=method, bounds_error=False
    )
    validation = rgi((gridlat, gridlon))
    # assert that the interpolated values are close to the expected values
    eps = np.finfo(np.float64).eps
    assert np.allclose(test, validation, atol=eps)
