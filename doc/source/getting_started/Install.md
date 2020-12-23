Installation
============

Presently geoid-toolkit is only available for use as a [GitHub repository](https://github.com/tsutterley/geoid-toolkit).
The contents of the repository can be download as a [zipped file](https://github.com/tsutterley/geoid-toolkit/archive/main.zip)  or cloned.
To use this repository, please fork into your own account and then clone onto your system.
```bash
git clone https://github.com/tsutterley/geoid-toolkit.git
```
Can then install using `setuptools`
```bash
python3 setup.py install
```
or `pip`
```bash
python3 -m pip install --user .
```

Alternatively can install the gravity_toolkit utilities from GitHub with `pip`
```bash
python3 -m pip install --user git+https://github.com/tsutterley/geoid-toolkit.git
```

Executable versions of this repository can also be tested using [Binder](https://mybinder.org/v2/gh/tsutterley/geoid-toolkit/main) and [Pangeo](https://binder.pangeo.io/v2/gh/tsutterley/geoid-toolkit/main).
