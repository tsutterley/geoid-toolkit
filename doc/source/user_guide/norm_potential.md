norm_potential.py
=================

- Calculates the normal potential at a given latitude and height
- Method described in [Barthelmes (2009)](http://icgem.gfz-potsdam.de/str-0902-revised.pdf) and [Moazezi (2012)](https://doi.org/10.1007/s12145-012-0102-2)

#### Calling Sequence
```python
from geoid_toolkit.norm_potential import norm_potential
U, dU_dr, dU_dtheta = norm_potential(lat, lon, h, 'WGS84', lmax)
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/norm_potential.py)

#### Inputs
1. `latitude`: latitude in degrees
2. `longitude`: longitude in degrees
3. `height`: height above reference ellipsoid in meters
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
5. `lmax`: maximum spherical harmonic degree

#### Outputs
- `U`: normal potential at height h
- `dU_dr`: derivative of normal potential with respect to radius
- `dU_dtheta`: derivative of normal potential with respect to theta