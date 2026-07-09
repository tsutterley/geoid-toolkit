# geoid-toolkit

Python tools for obtaining and working with static gravity field coefficients and calculating geoid heights

- [GFZ International Centre for Global Earth Models (ICGEM)](http://icgem.gfz-potsdam.de)

## About

<table>
  <tr>
    <td><b>Version:</b></td>
    <td>
        <a href="https://pypi.python.org/pypi/geoid-toolkit/" alt="PyPI"><img src="https://img.shields.io/pypi/v/geoid-toolkit.svg"></a>
        <a href="https://anaconda.org/conda-forge/geoid-toolkit" alt="conda-forge"><img src="https://img.shields.io/conda/vn/conda-forge/geoid-toolkit"></a>
        <a href="https://github.com/tsutterley/geoid-toolkit/releases/latest" alt="commits-since"><img src="https://img.shields.io/github/commits-since/tsutterley/geoid-toolkit/latest"></a>
    </td>
  </tr>
  <tr>
    <td><b>Citation:</b></td>
    <td>
        <a href="https://doi.org/10.5281/zenodo.5768711" alt="zenodo"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5768711.svg"></a>
    </td>
  </tr>
  <tr>
    <td><b>Tests:</b></td>
    <td>
        <a href="https://geoid-toolkit.readthedocs.io/en/latest/?badge=latest" alt="Documentation Status"><img src="https://readthedocs.org/projects/geoid-toolkit/badge/?version=latest"></a>
        <a href="https://github.com/tsutterley/geoid-toolkit/actions/workflows/ruff-format.yml" alt="Ruff"><img src="https://github.com/tsutterley/geoid-toolkit/actions/workflows/ruff-format.yml/badge.svg"></a>
    </td>
  </tr>
  <tr>
    <td><b>License:</b></td>
    <td>
        <a href="https://github.com/tsutterley/geoid-toolkit/blob/main/LICENSE" alt="License"><img src="https://img.shields.io/github/license/tsutterley/geoid-toolkit"></a>
    </td>
  </tr>
</table>

For more information: see the documentation at [geoid-toolkit.readthedocs.io](https://geoid-toolkit.readthedocs.io/)

## Installation

From PyPI:

```bash
python3 -m pip install geoid-toolkit
```

To include all optional dependencies:

```bash
python3 -m pip install geoid-toolkit[all]
```

Using `conda` or `mamba` from conda-forge:

```bash
conda install -c conda-forge geoid-toolkit
```

```bash
mamba install -c conda-forge geoid-toolkit
```

Development version from GitHub:

```bash
python3 -m pip install git+https://github.com/tsutterley/geoid-toolkit.git
```

### Running with Pixi

Alternatively, you can use [Pixi](https://pixi.sh/) for a streamlined workspace environment:

1. Install Pixi following the [installation instructions](https://pixi.sh/latest/#installation)
2. Clone the project repository:

```bash
git clone https://github.com/tsutterley/geoid-toolkit.git
```

3. Move into the `geoid-toolkit` directory

```bash
cd geoid-toolkit
```

4. Install dependencies and start JupyterLab:

```bash
pixi run start
```

This will automatically create the environment, install all dependencies, and launch JupyterLab in the [notebooks](./doc/source/notebooks/) directory.

## Resources

- [NASA GRACE mission site](https://www.nasa.gov/mission_pages/Grace/index.html)
- [NASA GRACE-FO mission site](https://www.nasa.gov/missions/grace-fo)
- [JPL GRACE Tellus site](https://grace.jpl.nasa.gov/)
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)

## Dependencies

- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)
- [numpy: Scientific Computing Tools For Python](https://www.numpy.org)
- [platformdirs: Python module for determining platform-specific directories](https://pypi.org/project/platformdirs/)
- [python-dateutil: powerful extensions to datetime](https://dateutil.readthedocs.io/en/stable/)

## Download

The program homepage is:  
<https://github.com/tsutterley/geoid-toolkit>

A zip archive of the latest version is available directly at:  
<https://github.com/tsutterley/geoid-toolkit/archive/main.zip>

## Disclaimer

This package includes software developed at the University of California at Irvine (UCI), the NASA Jet Propulsion Laboratory (JPL), NASA Goddard Space Flight Center (GSFC) and the University of Washington Applied Physics Laboratory (UW-APL).
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.
The software is provided here for your convenience but *with no guarantees whatsoever*.

## Contributing

This project contains work and contributions from the [scientific community](./CONTRIBUTORS.md).
If you would like to contribute to the project, please have a look at the [contribution guidelines](./doc/source/getting_started/Contributing.rst), [open issues](https://github.com/tsutterley/geoid-toolkit/issues) and [discussions board](https://github.com/tsutterley/geoid-toolkit/discussions).

## References

> Drewes, Kuglitsch, Ad&aacute;m and R&oacute;zsa
> "The Geodesist's Handbook 2016",
> Journal of Geodesy, 90, 907-1205 (2016).
> [doi: 10.1007/s00190-016-0948-z](https://doi.org/10.1007/s00190-016-0948-z)
>
> Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005).
> [doi: 10.1007/978-3-211-33545-1](https://doi.org/10.1007/978-3-211-33545-1)
>
> Holmes and Featherstone,
> "A Unified Approach to the Clenshaw Summation and the Recursive Computation of Very High Degree and Order Normalised Associated Legendre Functions",
> Journal of Geodesy (2002).
> [doi: 10.1007/s00190-002-0216-2](https://doi.org/10.1007/s00190-002-0216-2)
>
> Ince, Barthelmes, Rei&szlig;land, Elger, F&ouml;rste, Flechtner, and Schuh,
> "ICGEM -- 15 years of successful collection and distribution of global gravitational models, associated services, and future plans"
> Earth System Science Data, 11, 647--674 (2019).
> [doi: 10.5194/essd-11-647-2019](https://doi.org/10.5194/essd-11-647-2019)
>
> Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", (1981).
>
> Moazezi and Zomorrodian,
> "GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly",
> Earth Science Informatics (2012).
> [doi: 10.1007/s12145-012-0102-2](https://doi.org/10.1007/s12145-012-0102-2)
>
> Tscherning and Poder,
> "Some Geodetic Applications of Clenshaw Summation",
> Bollettino di Geodesia e Scienze, (1982).
>
> Wahr, Molenaar and Frank,
> "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE",
> Journal of Geophysical Research: Solid Earth, 103(B12), 30205-30229, (1998).
>  [doi: 10.1029/98JB02844](https://doi.org/10.1029/98JB02844)

## License

The content of this project is licensed under the [Creative Commons Attribution 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/) and the source code is licensed under the [MIT license](LICENSE).
