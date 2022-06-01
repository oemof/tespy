Thermal Engineering Systems in Python
=====================================
TESPy stands for "Thermal Engineering Systems in Python" and provides a
powerful simulation toolkit for thermal engineering plants such as power
plants, district heating systems or heat pumps. It is an external extension
module within the `Open Energy Modelling Framework <https://oemof.org/>`_ and
can be used as a standalone package.

.. figure:: https://raw.githubusercontent.com/oemof/tespy/9915f013c40fe418947a6e4c1fd0cd0eba45893c/docs/api/_images/logo_tespy_big.svg
    :align: center

With the TESPy package you are able to calculate stationary operation in order
to design the process of thermal energy systems. From that point it is possible
to simulate the offdesign behavior of your plant using underlying
characteristics for each of the plants components. The package includes basic
components, such as turbines, pumps, compressors, heat exchangers, pipes,
mixers and splitters as well as some advanced components (derivatives of heat
exchangers, drum).

Everybody is welcome to use and/or develop TESPy. Contribution is already
possible on a low level by simply fixing typos in TESPy's documentation or
rephrasing sections which are unclear. If you want to support us that way
please fork the TESPy repository to your own github account and make changes
as described in the github guidelines:
https://guides.github.com/activities/hello-world/

Key Features
============
* **Open** Source
* **Generic** thermal engineering applications
* **Automatic** model documentation in LaTeX for high transparency and
  reproducibility
* **Extendable** framework for the implementation of custom components and
  component groups
* **Postprocessing** features like exergy analysis and fluid property plotting

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - |pytests| |checks| |packaging| |coveralls|
    * - package
      - | |version| |wheel| |supported-versions| |commits-since|
    * - reference
      - |joss| |zenodo|

.. |docs| image:: https://readthedocs.org/projects/tespy/badge/?style=flat
    :target: https://readthedocs.org/projects/tespy
    :alt: Documentation Status

.. |pytests| image:: https://github.com/oemof/tespy/workflows/tox%20pytests/badge.svg
    :target: https://github.com/oemof/tespy/actions?query=workflow%3A%22tox+pytests%22
    :alt: tox pytest

.. |checks| image:: https://github.com/oemof/tespy/workflows/tox%20checks/badge.svg
    :target: https://github.com/oemof/tespy/actions?query=workflow%3A%22tox+checks%22
    :alt: tox checks

.. |packaging| image:: https://github.com/oemof/tespy/workflows/packaging/badge.svg
    :target: https://github.com/oemof/tespy/actions?query=workflow%3Apackaging
    :alt: packaging

.. |coveralls| image:: https://coveralls.io/repos/oemof/tespy/badge.svg?branch=main&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/oemof/tespy

.. |version| image:: https://img.shields.io/pypi/v/tespy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/tespy

.. |wheel| image:: https://img.shields.io/pypi/wheel/tespy.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/tespy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/tespy.svg
    :alt: Supported Python versions
    :target: https://pypi.org/project/tespy

.. |commits-since| image:: https://img.shields.io/github/commits-since/oemof/tespy/latest/dev
    :alt: Commits since latest release
    :target: https://github.com/oemof/tespy/compare/main...dev

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2555866.svg
   :alt: Release archive
   :target: https://doi.org/10.5281/zenodo.2555866

.. |joss| image:: https://joss.theoj.org/papers/590b0b4767606bce4d0ebe397d4b7a4f/status.svg
   :alt: Software Paper in JOSS
   :target: https://joss.theoj.org/papers/590b0b4767606bce4d0ebe397d4b7a4f

.. end-badges

Documentation
=============
You can find the full documentation at
`readthedocs <http://tespy.readthedocs.org>`_. Use the
`project site <http://readthedocs.org/projects/tespy>`_ of readthedocs to
choose the version of the documentation. Go to the
`download page <http://readthedocs.org/projects/tespy/downloads/>`_ to
download different versions and formats (pdf, html, epub) of the documentation.

To get the latest news visit and follow our `website <https://www.oemof.org>`_.

Installing TESPy
================
If you have a working Python3 environment, use pypi to install the latest
tespy version:

.. code:: bash

  pip install tespy

If you want to use the latest features, you might want to install the
**developer version**. See section
`Developing TESPy <http://tespy.readthedocs.io/en/dev/developing_tespy.html>`_
for more information. The developer version is not recommended for productive
use.

Get in touch
============

Online meetings
---------------
We have decided to start a reoccurring meeting for all interested TESPy users
and (potential) developers. You are invited to join us on every 3rd Monday of a
month at 17:00 CE(S)T. The first meeting will be held at June, 20, 2022.

The intent of this meeting is to establish a more active and well connected
network of TESPy users and developers. Participants will have a space to present
and discuss their models, get help and exchange experience with others, learn
about new features or suggest improvements or new features.

The preliminary schedule is as follows:

* 17:00 – 17:05: Short intro
* 17:05 – 17:50: User presentations and discussion
* 17:50 – 18:00: Break
* 18:00 – 18:40: New feature presentation/discussion, tutorial session
* 18:40 – 18:45: Feedback

.. note::

    The meetings will thrive through active participation. If you have
    suggestions for topics - **from beginner tutorial to suggestion or**
    **discussion of specific new features** - please let us know your interest
    when registering to the meeting. We appreciate your contribution!

To get a better overview on the number of participants and to plan the agenda
of each meeting accordingly, **please send an e-mail to**
**francesco.witte(at)dlr.de for registration**.

You can download the meeting series here:
:download:`regular_meeting.ics </regular_meeting.ics>`.

We are looking forward to meeting you!

User forum
----------
We have implemented a
`discussion room on GitHub <https://github.com/oemof/tespy/discussions>`__ as
user forum. If you have issues with setting up your model or any other question
about using the software, you are invited to start a discussion there.

Examples
========

For a short introduction on how TESPy works and how you can use it, we provide
some
`examples and tutorials <https://tespy.readthedocs.io/en/main/tutorials_examples.html>`_,
go and check them out. You can download the python scripts of all example plants
from the
`tespy_examples <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy>`_
repository.

Citation
========
The scope and functionalities of TESPy have been documented in a paper
published in the Journal of Open Source Software with an OpenAccess license.
Download the paper from https://doi.org/10.21105/joss.02178. As TESPy is a free
software, we kindly ask that you add a reference to TESPy if you use the
software for your scientific work. Please cite the article with the BibTeX
citation below.

BibTeX citation::

    @article{Witte2020,
        doi = {10.21105/joss.02178},
        year = {2020},
        publisher = {The Open Journal},
        volume = {5},
        number = {49},
        pages = {2178},
        author = {Francesco Witte and Ilja Tuschy},
        title = {{TESPy}: {T}hermal {E}ngineering {S}ystems in {P}ython},
        journal = {Journal of Open Source Software}
    }

Furthermore, a paper on the exergy analysis feature has been published in
the mdpi journal energies. You can download the pdf at
https://doi.org/10.3390/en15114087. If you are using this feature specifically,
you can reference it with the following BibTeX citation:

BibTeX citation::

    @article{Witte2022,
        doi = {10.3390/en15114087},
        year = {2022},
        publisher = {The Open Journal},
        volume = {15},
        number = {11},
        article-number = {4087},
        issn = {1996-1073},
        author = {Witte, Francesco and Hofmann, Mathias and Meier, Julius and Tuschy, Ilja and Tsatsaronis, George},
        title = {Generic and Open-Source Exergy Analysis&mdash;Extending the Simulation Framework TESPy},
        journal = {Energies}
    }


Additionally, you have the possibility to cite a specific version of TESPy to
make your work reproducible. The source code of every version is published on
zenodo. Find your version here: https://doi.org/10.5281/zenodo.2555866.

License
=======
Copyright (c) 2017-2022 oemof developer group

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
