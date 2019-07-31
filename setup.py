from setuptools import setup, find_packages
setup(
	name='geoid-toolkit',
	version='1.0.0.0',
	description='Reads gravity model coefficients and calculates geoid heights',
	url='https://github.com/tsutterley/geoid-toolkit',
	author='Tyler Sutterley',
	author_email='tsutterl@uw.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Physics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	keywords='static gravity field, geoid height',
	packages=find_packages(),
	install_requires=['numpy'],
)
