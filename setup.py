from setuptools import setup, find_packages

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# get install requirements
with open('requirements.txt') as fh:
    install_requires = fh.read().splitlines()

setup(
    name='geoid-toolkit',
    version='1.0.0.0',
    description='Reads gravity model coefficients and calculates geoid heights',
    long_description=long_description,
    long_description_content_type="text/markdown",
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
        'Programming Language :: Python :: 3.7',
    ],
    keywords='static gravity field, geoid height',
    packages=find_packages(),
    install_requires=install_requires,
    include_package_data=True,
)
