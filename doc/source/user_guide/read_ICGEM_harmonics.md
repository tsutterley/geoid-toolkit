read_ICGEM_harmonics.py
=======================

- Reads gravity model files and extracts spherical harmonic data from the [GFZ International Centre for Global Earth Models (ICGEM)](http://icgem.gfz-potsdam.de/)

#### Calling Sequence
```python
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
Ylms = read_ICGEM_harmonics(model_file)
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/read_ICGEM_harmonics.py)

#### Arguments
1. `model_file`: full path to GFZ ICGEM gfc spherical harmonic data file

#### Keyword arguments
- `LMAX`: maximum degree and order of output spherical harmonic coefficients
- `TIDE`: [tide system of output geoid](http://mitgcm.org/~mlosch/geoidcookbook/node9.html)
    * `'tide_free'`: no permanent direct and indirect tidal potentials
    * `'mean_tide'`: restores permanent tidal potentials (direct and indirect)
    * `'zero_tide'`: restores permanent direct tidal potential
- `FLAG`: string denoting data lines

#### Returns
- `clm`: cosine spherical harmonics of input data
- `slm`: sine spherical harmonics of input data
- `eclm`: cosine spherical harmonic standard deviations of type errors
- `eslm`: sine spherical harmonic standard deviations of type errors
- `modelname`: name of the gravity model
- `earth_gravity_constant`: GM constant of the Earth for the gravity model
- `radius`: semi-major axis of the Earth for the gravity model
- `max_degree`: maximum degree and order for the gravity model
- `errors`: error type of the gravity model
- `norm`: normalization of the spherical harmonics
- `tide_system`: tide system of gravity model (mean_tide, zero_tide, tide_free)