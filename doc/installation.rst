.. _installation_and_setup_label:

######################
Installation and setup
######################

.. contents::
    :depth: 1
    :local:
    :backlinks: top


Following you find guidelines for the installation process for different
operation systems.

Linux
=====

Installing Python 3
-------------------

TESPy is a Python package, thus it requires you to have Python 3 installed.
Most Linux distributions will have Python 3 in their repository. Use the
specific software management to install it. If you are using Ubuntu/Debian try
executing the following code in your terminal:

.. code:: console

  sudo apt-get install python3

You can also download different versions of Python via
https://www.python.org/downloads/.

Having Python 3 installed
-------------------------

.. code:: console

  pip install tespy

To use pip you have to install the pypi package. Normally pypi is part of your
virtual environment.

.. _virtualenv_label:

Using virtualenv instead of system wide Python
----------------------------------------------

Instead of installing TESPy with pip to your system Python, you can instead
install TESPy to a virtual Python environment.

 1. Install virtualenv using the package management of your Linux distribution,
    pip install or install it from source
    (`see virtualenv documentation <https://virtualenv.pypa.io/en/stable/installation/>`_)
 2. Open terminal to create and activate a virtual environment by typing:

    .. code-block:: console

       virtualenv -p /usr/bin/python3 your_env_name
       source your_env_name/bin/activate

 3. In terminal type: :code:`pip install tespy`

Warning: If you have an older version of virtualenv you should update pip
:code:`pip install --upgrade pip`.

.. _tespy_installation_windows_label:

Windows
=======

Having Python 3 installed
-------------------------

If you already have a working Python 3 environment you can install TESPy by
using pip. We recommend you installing the package in a virtual environment.
You can use virtualenv (:ref:`see here for instructions <virtualenv_label>`)
or a virtual environment e. g. in :ref:`Anaconda <anaconda_label>`.

TESPy requires the CoolProp python package, which will be installed
automatically. As it is build from C++, you need to have Microsoft Visual
Studio C++ Build Tools 14.0 installed. Try installing TESPy with the following
code in the command window of your python environment:

.. code:: console

  pip install tespy

If you get an error within the installation of the CoolProp package, install
the C++ Build Tools first and then restart the installation. Also, if pip is
not part of your python environment, you have to install the pypi package.

.. _anaconda_label:

Using Anaconda
--------------

 1. Download latest `Anaconda <https://www.continuum.io/downloads#windows>`_
    for Python 3.x (64 or 32 bit).
 2. Install Anaconda
 3. Open 'Anaconda Prompt' to create and activate a virtual environment by
    typing:

    .. code-block:: console

       conda create -n yourenvname python=3.x
       activate yourenvname

 4. In the active Anaconda Prompt type: :code:`pip install tespy`
 5. If the installation of CoolProp fails, make shure, you have
    Microsoft Visual Stuido C++ Build Tools 14.0 installed on your computer.


Mac OSX
=======

Installation instructions for Mac OSX are not available, yet. If you want to
share your knwolegde on the installation and fill this gap, feel free to
contact us.
