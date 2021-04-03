import os
from setuptools import setup, find_packages

# get long_description from README.rst
with open("README.rst", "r") as fh:
    long_description = fh.read()

# get install requirements
with open('requirements.txt') as fh:
    install_requires = fh.read().splitlines()

# get version
with open('version.txt') as fh:
    version = fh.read()

# list of all scripts to be included with package
scripts=[os.path.join('scripts',f) for f in os.listdir('scripts') if f.endswith('.py')]

setup(
    name='geoid-toolkit',
    version=version,
    description='Reads gravity model coefficients and calculates geoid heights',
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url='https://github.com/tsutterley/geoid-toolkit',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='static gravity field, geoid height',
    packages=find_packages(),
    install_requires=install_requires,
    scripts=scripts,
    include_package_data=True,
)
