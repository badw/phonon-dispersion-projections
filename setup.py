
"""
phonon-dispersion-projections
"""

from os.path import abspath, dirname
from setuptools import find_packages, setup

setup(
    name='pdp',
    version='1.0.0',
    description='plotting phonon dispersion projections',
    url="https://github.com/badw/phonon-dispersion-projections",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    license='MIT',
    packages=find_packages(),
    install_requires=['pymatgen','numpy','sumo','phonopy','matplotlib']
    )
