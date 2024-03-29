#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.md').read()
doclink = """
Documentation
-------------

The full documentation is at http://dreem_tools.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='dreem_tools',
    version='0.1.0',
    description='a set of tools for setting up and processing data generated from dreem',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/dreem_tools',
    packages=[
        'dreem_tools',
    ],
    package_dir={'dreem_tools': 'dreem_tools'},
    py_modules=[
        'dreem_tools/dataframe',
        'dreem_tools/cli',
        'dreem_tools/logger',
        'dreem_tools/motif',
        'dreem_tools/plotting',
        'dreem_tools/process',
        'dreem_tools/parse',
        'dreem_tools/run',
        'dreem_tools/util'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='dreem_tools',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
            'dreem-tools = dreem_tools.cli : cli',
        ]
    }
)
