real_potential.py
=================

- Calculates the real potential at a given latitude and height using coefficients from a gravity model
- Method described in [Barthelmes (2009)](http://icgem.gfz-potsdam.de/str-0902-revised.pdf) and [Moazezi (2012)](https://doi.org/10.1007/s12145-012-0102-2)

#### Calling Sequence
```python
from gravity_toolkit.real_potential import real_potential
W, dW_dr, dW_dtheta = real_potential(lat, lon, h, clm, slm, lmax)
```
[Source code](https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/real_potential.py)

#### Inputs
1. `latitude`: latitude in degrees
2. `longitude`: longitude in degrees
3. `height`: height above reference ellipsoid in meters
4. `GM`: geocentric graviational constant used in gravity model
5. `R`: average radius used in gravity model
6. `clm`: cosine spherical harmonics for a gravity model
7. `slm`: sine spherical harmonics for a gravity model
8. `lmax`: maximum spherical harmonic degree

#### Options
- `GAUSS`: Gaussian Smoothing Radius in km

#### Outputs
- `W`: real potential at height h
- `dW_dr`: derivative of real potential with respect to radius
- `dW_dtheta`: derivative of real potential with respect to theta