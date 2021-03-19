norm_gravity.py
===============

- Calculates the normal gravity of an ellipsoid at a given latitude and height and calculates the derivative with respect to height

#### Calling Sequence
```python
from geoid_toolkit.norm_gravity import norm_gravity
gamma_h,dgamma_dh = norm_gravity(latitude, height, 'WGS84')
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/norm_gravity.py)

#### Arguments
1. `latitude`: latitude in degrees
2. `height`: height above reference ellipsoid in meters
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

#### Returns
- `gamma_h`: normal gravity for ellipsoid at height
- `dgamma_dh`: derivative of normal gravity with respect to height