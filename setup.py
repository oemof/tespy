#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='TESPy',
    version='0.6.2',
    license='MIT',
    description='Thermal Engineering Systems in Python (TESPy)',
    long_description='%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub(
            '', read('README.rst')
        )
    ),
    author='Francesco Witte',
    author_email='francesco.witte@dlr.de',
    url='https://github.com/oemof/tespy',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    package_data={'src': ['*.json']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering',
    ],
    project_urls={
        'Documentation': 'https://tespy.readthedocs.io/',
        'Changelog': 'https://tespy.readthedocs.io/en/main/whats_new.html',
        'Issue Tracker': 'https://github.com/oemof/tespy/issues',
    },
    python_requires='>=3.7',
    install_requires=[
        'CoolProp>=6.4,<7',
        'matplotlib>=3.2.1,<4',
        'numpy>=1.13.3,<2',
        'pandas>=1.3.0,<2',
        'tabulate>=0.8.2,<0.9'
    ],
    extras_require={
        'dev': [
            'furo',
            'py',
            'pytest',
            'sphinx',
            'sphinx-copybutton',
            'sphinx-design',
            'sphinxcontrib.bibtex',
            'tox',
        ],
        'dummy': ['tespy']
    }
)
