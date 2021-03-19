gravity_anomaly.py
==================

- Calculates the gravity anomaly at a given latitude and longitude using different methods
- Method described in [Barthelmes (2009)](http://icgem.gfz-potsdam.de/str-0902-revised.pdf), [Moazezi (2012)](https://doi.org/10.1007/s12145-012-0102-2) and [Molodensky (1958)](https://doi.org/10.1007/BF02537957)

#### Calling Sequence
```python
from geoid_toolkit.gravity_anomaly import gravity_anomaly
zeta = gravity_anomaly(lat, lon, h, 'WGS84', clm, slm, lmax, R, GM, METHOD='first')
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/height_anomaly.py)

#### Inputs
1. `lat`: latitudinal points to calculate geoid height
2. `lon`: longitudinal points to calculate geoid height
3. `h`: ellipsoidal height
4. `refell`: reference ellipsoid name
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
5. `clm`: cosine spherical harmonics for a gravity model
6. `slm`: sine spherical harmonics for a gravity model
7. `lmax`: maximum spherical harmonic degree
8. `R`: average radius used in gravity model
9. `GM`: geocentric graviational constant used in gravity model

#### Options
- `METHOD`: method for calculating gravity anomalies
    * `first`: classic first approximation method
    * `second`: classic second approximation method
    * `molodensky`: Molodensky method (1958)
- `GAUSS`: Gaussian Smoothing Radius in km

#### Outputs
- `ddelta_g`: gravity anomaly for a given ellipsoid in meters
