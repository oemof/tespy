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
    name='tespy',
    version='0.3.0 dev',
    license='MIT',
    description='Thermal Engineering Systems in Python',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub(
            '', read('README.rst')
        ),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Francesco Witte',
    author_email='francesco.witte@hs-flensburg.de',
    url='https://github.com/oemof/tespy',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Utilities',
    ],
    project_urls={
        'Documentation': 'https://tespy.readthedocs.io/',
        'Changelog': 'https://tespy.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/oemof/tespy/issues',
    },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.6, <3.8',
    install_requires=[
        'CoolProp>=6,<7',
        'numpy>=1.13.3,<2',
        'pandas>=0.19.2,!=1.0.0,<2',
        'scipy>=0.19.1,<2',
        'tabulate>=0.8.2,<0.9'
    ],
    extras_require={'dev': ['pytest', 'sphinx', 'sphinx_rtd_theme', ],
                    'dummy': ['tespy']}
)
