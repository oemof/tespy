.. _installation_and_setup_label:

######################
Installation and setup
######################

.. contents::
    :depth: 1
    :local:
    :backlinks: top


Following you find guidelines for the installation process for different operation systems.

Linux
=====

Having Python 3 installed
-------------------------

TESPy is a Python package, thus it requires you to have Python 3 installed. It is highly recommended to
use a virtual environment. See this `tutorial <https://docs.python.org/3/tutorial/venv.html>`_ for more
help or see the sections below. If you already have a Python 3 environment you can clone the
`github repository <https://github.com/oemof/tespy>`_ and install it:

.. code:: console

  pip install -e /path/to/the/repository

To use pip you have to install the pypi package. Normally pypi is part of your virtual environment.

Using Linux repositories to install Python
------------------------------------------

Most Linux distributions will have Python 3 in their repository. Use the specific software management to install it. 
If you are using Ubuntu/Debian try executing the following code in your terminal: 

.. code:: console

  sudo apt-get install python3
  
You can also download different versions of Python via https://www.python.org/downloads/.

Using Virtualenv (community driven)
-----------------------------------

Skip the steps you have already done. Check your architecture first (32/64 bit).

 1. Install virtualenv using the package management of your Linux distribution, pip install or install
    it from source (`see virtualenv documentation <https://virtualenv.pypa.io/en/stable/installation/>`_)
 2. Open terminal to create and activate a virtual environment by typing:

    .. code-block:: console

       virtualenv -p /usr/bin/python3 your_env_name
       source your_env_name/bin/activate

 3. In terminal type: :code:`pip install -e /path/to/the/repository`
 
Warning: If you have an older version of virtualenv you should update pip :code:`pip install --upgrade pip`.

Windows
=======

Installation instructions for Windows are not available, yet. If you want to share your knwolegde on the installation and fill this gap, feel free to contact us.

Mac OSX
=======

Installation instructions for Mac OSX are not available, yet. If you want to share your knwolegde on the installation and fill this gap, feel free to contact us.