.. image:: https://coveralls.io/repos/github/oemof/tespy/badge.svg?branch=dev
    :target: https://coveralls.io/github/oemof/tespy?branch=dev
.. image:: https://travis-ci.org/oemof/tespy.svg?branch=dev
    :target: https://travis-ci.org/oemof/tespy
.. image:: https://readthedocs.org/projects/tespy/badge/?version=dev
    :target: https://tespy.readthedocs.io/en/dev/
.. image:: https://badge.fury.io/py/TESPy.svg
    :target: https://badge.fury.io/py/TESPy
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2555866.svg
   :target: https://doi.org/10.5281/zenodo.2555866

TESPy stands for "Thermal Engineering Systems in Python" and provides a powerful simulation toolkit for thermal engineering plants such as power plants, district heating systems or heat pumps. It is an external extension module within the `Open Energy Modelling Framework <https://oemof.org/>`_ and can be used as a standalone package.

With the TESPy package you are able to calculate stationary operation in order to design the process of thermal energy systems. From that point it is possible to simulate the offdesign behaviour of your plant using underlying characteristics for each of the plants components. The package inlcudes basic components, such as turbines, pumps, compressors, heat exchangers, pipes, mixers and splitters as well as some advanced components (derivatives of heat exchangers, drum).

Everybody is welcome to use and/or develop TESPy. Contribution is already possible on a low level by simply fixing typos in TESPy's documentation or rephrasing sections which are unclear. If you want to support us that way please fork the TESPy repository to your own github account and make changes as described in the github guidelines: https://guides.github.com/activities/hello-world/

Documentation
=============

You can find the full documentation at `readthedocs <http://tespy.readthedocs.org>`_. Use the `project site <http://readthedocs.org/projects/tespy>`_ of readthedocs to choose the version of the documentation. Go to the `download page <http://readthedocs.org/projects/tespy/downloads/>`_ to download different versions and formats (pdf, html, epub) of the documentation.

To get the latest news visit and follow our `website <https://www.oemof.org>`_.

Installing TESPy
================

If you have a working Python3 environment, use pypi to install the latest tespy version:

.. code:: bash

  pip install tespy

If you want to use the latest features, you might want to install the **developer version**. See section `Developing TESPy <http://tespy.readthedocs.io/en/latest/developing_tespy.html>`_ for more information. The developer version is not recommended for productive use.

Examples
========

For a short introduction on how TESPy works and how you can use it, we provide some `basic examples and tutorials <http://tespy.readthedocs.io/en/latest/getting_started.html>`_, go and check them out. You can download the python scripts of all example plants from the `tespy_examples repository <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy>`_.

License
=======

Copyright (c) 2019 oemof developer group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
