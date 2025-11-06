~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Thermal Engineering Systems in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TESPy stands for "Thermal Engineering Systems in Python" and provides a
powerful simulation toolkit for thermal engineering plants such as various
types of power plants (including organic rankine cycles), heat pumps or
refrigeration machines. Due to its flexibility it is actually possible to
model any kind of thermal energy conversion process, this also includes energy
balancing of industrial processes, district heating or HVAC systems. It is
part of the Open Energy Modelling Framework `oemof <https://oemof.org/>`_ and
can be used as a standalone package.

.. image:: /_static/images/logo_tespy_big.svg
   :align: center
   :class: only-light


.. image:: /_static/images/logo_tespy_big_darkmode.svg
   :align: center
   :class: only-dark

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
please fork the `TESPy repository <https://github.com/oemof/tespy>`_ to your
own GitHub account and make changes as described in the GitHub guidelines:
https://guides.github.com/activities/hello-world/

Key Features
============
* **Free** and **Open** Source Software
* **Flexible** models of thermal engineering applications due to component
  based architecture
* **Extendable** framework for the implementation of custom components, fluid
  property formulations and equations
* **Integration** of optimization capabilities through an API to pygmo
* **Postprocessing** features like exergy analysis and fluid property plotting

Quick installation
==================

If you have a working Python3 environment, use pypi to install the latest
tespy version.

.. code:: bash

    pip install tespy

We provide more detailed
:ref:`installation instructions <installation_and_setup_label>`, too.

If you want to use the latest features, you might want to install the
**developer version**. See
:ref:`this section <tespy_development_how_label>` for more information.

Getting into TESPy
==================

For a good start on how TESPy works and how you can use it, we provide some
:ref:`basic <tespy_basics_label>` and :ref:`advanced <tespy_tutorial_label>`
tutorials in the User Guide section. The
:ref:`modules <tespy_modules_label>` section provides you with in depth
information on the different modules of TESPy.

Citation
========

The scope and functionalities of TESPy have been documented in a paper
published in the Journal of Open Source Software with an Open-Access license.
Download the paper from https://doi.org/10.21105/joss.02178 :cite:`Witte2020`.
As TESPy is a free software, we kindly ask that you add a reference to TESPy
if you use the software for your scientific work. Please cite the article with
the BibTeX citation below.

BibTeX citation

.. code::

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

Additionally, you have the possibility to cite a specific version of TESPy to
make your work reproducible. The source code of every version is published on
zenodo. Find your version here: https://doi.org/10.5281/zenodo.2555866.

License
=======

.. include:: ../LICENSE
