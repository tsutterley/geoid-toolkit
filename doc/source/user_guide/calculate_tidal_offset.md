calculate_tidal_offset.py
=========================

- Calculates the spherical harmonic offset for a tide system to change from a tide free state where there is no permanent direct and indirect tidal potentials

#### Calling Sequence
```python
from geoid_toolkit.calculate_tidal_offset import calculate_tidal_offset
delta = calculate_tidal_offset(TIDE, GM, R, refell)
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/calculate_tidal_offset.py)

#### Arguments
1. `TIDE`: output tidal system
    `'mean_tide'`: restores permanent tidal potentials (direct and indirect)
    `'zero_tide'`: restores permanent direct tidal potential
2. `R`: average radius used in gravity model
3. `GM`: geocentric graviational constant used in gravity model
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

#### Returns
- `delta`: offset for changing from tide free system
