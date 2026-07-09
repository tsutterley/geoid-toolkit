import os
from setuptools import setup

# list of all scripts to be included with package
scripts = [
    os.path.join('scripts', f)
    for f in os.listdir('scripts')
    if f.endswith('.py')
]

setup(
    name='geoid-toolkit',
    scripts=scripts,
)
