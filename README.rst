========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor|
        | |coveralls| |codecov|
        | |codacy| |codeclimate|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/tespy/badge/?style=flat
    :target: https://readthedocs.org/projects/tespy
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/oemof/tespy.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/oemof/tespy

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/oemof/tespy?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/oemof/tespy

.. |coveralls| image:: https://coveralls.io/repos/oemof/tespy/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/oemof/tespy

.. |codecov| image:: https://codecov.io/gh/oemof/tespy/branch/master/graphs/badge.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/oemof/tespy

.. |codacy| image:: https://img.shields.io/codacy/grade/CODACY_TOKEN.svg
    :target: https://www.codacy.com/app/oemof/tespy
    :alt: Codacy Code Quality Status

.. |codeclimate| image:: https://codeclimate.com/github/oemof/tespy/badges/gpa.svg
   :target: https://codeclimate.com/github/oemof/tespy
   :alt: CodeClimate Quality Status

.. |version| image:: https://img.shields.io/pypi/v/tespy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/tespy

.. |wheel| image:: https://img.shields.io/pypi/wheel/tespy.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/tespy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/tespy.svg
    :alt: Supported versions
    :target: https://pypi.org/project/tespy

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/tespy.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/tespy

.. |commits-since| image:: https://img.shields.io/github/commits-since/oemof/tespy/v0.3.0 dev.svg
    :alt: Commits since latest release
    :target: https://github.com/oemof/tespy/compare/v0.3.0 dev...master



.. end-badges

Thermal Engineering Systems in Python

* Free software: MIT license

Installation
============

::

    pip install tespy

You can also install the in-development version with::

    pip install https://github.com/oemof/tespy/archive/master.zip


Documentation
=============


https://tespy.readthedocs.io/


Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
