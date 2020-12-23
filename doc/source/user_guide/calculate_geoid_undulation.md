calculate_geoid_undulation.py
=============================

- Wrapper function for computing geoid undulations from a gravity model
- Method described in [Barthelmes (2009)](http://icgem.gfz-potsdam.de/str-0902-revised.pdf) and [Moazezi (2012)](https://doi.org/10.1007/s12145-012-0102-2)

#### Calling Sequence
```python
from gravity_toolkit.calculate_geoid_undulation import calculate_geoid_undulation
N = calculate_geoid_undulation(lon, lat, gravity_model_file, REFERENCE='WGS84')
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/calculate_geoid_undulation.py)

#### Inputs
1. `lon`: longitudinal points to calculate geoid height
2. `lat`: latitudinal points to calculate geoid height
3. `gravity_model_file`: full path to static gravity file with spherical harmonic coefficients

#### Options
- `REFERENCE`: reference ellipsoid name
    * `'CLK66'`: Clarke 1866
    * `'GRS67'`: Geodetic Reference System 1967
    * `'GRS80'`: Geodetic Reference System 1980
    * `'WGS72'`: World Geodetic System 1972
    * `'WGS84'`: World Geodetic System 1984
    * `'ATS77'`: Quasi-earth centred ellipsoid for ATS77
    * `'NAD27'`: North American Datum 1927 (=CLK66)
    * `'NAD83'`: North American Datum 1983 (=GRS80)
    * `'INTER'`: International
    * `'KRASS'`: Krassovsky (USSR)
    * `'MAIRY'`: Modified Airy (Ireland 1965/1975)
    * `'TOPEX'`: TOPEX/POSEIDON ellipsoid
    * `'EGM96'`: EGM 1996 gravity model
- `LMAX`: maximum spherical harmonic degree
- `TIDE`: tide system of output geoid
- `GAUSS`: Gaussian Smoothing Radius in km

#### Outputs
- `N`: geoidal undulation for a given ellipsoid in meters