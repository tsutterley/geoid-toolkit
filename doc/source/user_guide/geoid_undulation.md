geoid_undulation.py
===================

- Calculates the geoidal undulation at a given latitude and longitude using an iterative approach
- Method described in [Barthelmes (2009)](http://icgem.gfz-potsdam.de/str-0902-revised.pdf) and [Moazezi (2012)](https://doi.org/10.1007/s12145-012-0102-2)

#### Calling Sequence
```python
from geoid_toolkit.geoid_undulation import geoid_undulation
N = geoid_undulation(lat, lon, 'WGS84', clm, slm, lmax, R, GM)
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/geoid_undulation.py)

#### Arguments
1. `lat`: latitudinal points to calculate geoid height
2. `lon`: longitudinal points to calculate geoid height
3. `refell`: reference ellipsoid name
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
4. `clm`: cosine spherical harmonics for a gravity model
5. `slm`: sine spherical harmonics for a gravity model
6. `lmax`: maximum spherical harmonic degree
7. `R`: average radius used in gravity model
8. `GM`: geocentric graviational constant used in gravity model

#### Keyword arguments
- `GAUSS`: Gaussian Smoothing Radius in km
- `EPS`: level of precision for calculating geoid height

#### Returns
- `N`: geoidal undulation for a given ellipsoid in meters