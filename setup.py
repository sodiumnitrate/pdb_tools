from setuptools import setup, find_packages
import pdb_tools
import sys
import os

setup(
        name = 'pdb_tools',
        version = pdb_tools.__version__,
        description = 'some protein scripts/tools for fun',
        url = 'https://github.com/sodiumnitrate/pdb_tools.git',
        author = 'Irem Altan',
        author_email = 'irem.altan@yale.edu',
        packages = find_packages(),
        install_requires = ['gemmi','pypdb','numpy','scikit-spatial'],
        python_requires = '>3.6',
        )
