.. _tespy_label:

~~~~~~~~~~~
About TESPy
~~~~~~~~~~~

TESPy stands for "Thermal Engineering Systems in Python" and provides a
powerful simulation toolkit for thermal engineering plants such as power
plants, district heating systems or heat pumps. It is an external extension
module within the `Open Energy Modeling Framework <https://oemof.org/>`_ and
can be used as a standalone package.

With the TESPy package you are able to calculate stationary operation in order
to design the process of your plant. From that point it is possible to
predict the offdesign behavior of your plant using underlying characteristics
for each of the plants components. For now, the package includes basic
components, such as turbines, pumps, compressors, heat exchangers, pipes,
mixers and splitters as well as some advanced components
(derivatives of heat exchangers, drum).

Everybody is welcome to use and/or develop TESPy. Contribution is already
possible on a low level by simply fixing typos in TESPy's documentation or
rephrasing sections which are unclear. If you want to support us that way
please fork the TESPy repository to your own github account and make
changes as described in the github guidelines:
https://guides.github.com/activities/hello-world/

.. contents::
    :depth: 1
    :local:
    :backlinks: top

Quick installation
==================

If you have a working Python3 environment, use pypi to install the latest
tespy version.

.. code:: bash

    pip install tespy

We provide more detailed
:ref:`installation instructions <installation_and_setup_label>`, too.

If you want to use the latest features, you might want to install the
**developer version**. See section
:ref:`Developing TESPy <developing_tespy_label>` for more information.

Getting into TESPy
==================

For a good start on how TESPy works and how you can use it, we provide some
examples at the :ref:`examples section <tespy_examples_label>`. On top of that,
there is a step-by-step :ref:`tutorial <heat_pump_tutorial_label>` on how to
model a heat pump in TESPy.

The :ref:`TESPy modules section <tespy_modules_label>` provides you
with all important information on the different modules of TESPy.

Citation
========

The scope and functionalities of TESPy have been documented in a paper
published in the Journal of Open Source Software with an OpenAccess license.
Download the paper from https://doi.org/10.21105/joss.02178. As TESPy is a free
software, we kindly ask that you add a reference to TESPy if you use the
software for your scientific work. Please cite the article with the BibTeX
citation below.

Additionally, you have the possibility to cite a specific version of TESPy to
make your work reproducible. The source code of every version is published on
zenodo. Find your version here: https://doi.org/10.5281/zenodo.2555866.

BibTeX citation::

    @article{Witte2020,
        doi = {10.21105/joss.02178}
        year = {2020},
        publisher = {The Open Journal},
        volume = {5},
        number = {49},
        pages = {2179},
        author = {Francesco Witte and Ilja Tuschy},
        title = {{TESPy}: {T}hermal {E}ngineering {S}ystems in {P}ython},
        journal = {Journal of Open Source Software}
    }

License
=======

Copyright (c) 2017-2020 oemof developer group

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
