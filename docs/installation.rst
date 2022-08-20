.. _installation_and_setup_label:

######################
Installation and setup
######################

Following you find guidelines for the installation process for linux and
windows. TESPy is a Python package, thus it requires you to have Python 3
installed.

.. tab-set::

   .. tab-item:: Linux

      **Installing Python 3**

      Most Linux distributions will have Python 3 in their repository. Use the
      specific software management to install it, if it is not yet installed. If
      you are using Ubuntu/Debian try executing the following code in your
      terminal:

      .. code:: console

         sudo apt-get install python3

      You can also download different versions of Python via
      https://www.python.org/downloads/.

      **Having Python 3 installed**

      We recommend installting TESPy within a virtual Python enviroment an not
      into the base, system wide Python installation. On Linux you can use
      virtualenv to do so.

      1. Install virtualenv using the package management of your Linux distribution,
         pip install or install it from source
         (`see virtualenv documentation <https://virtualenv.pypa.io/en/stable/installation.html>`_)
      2. Open terminal to create and activate a virtual environment by typing:

         .. code-block:: console

            virtualenv -p /usr/bin/python3 your_env_name
            source your_env_name/bin/activate

      3. In terminal type: :code:`pip install tespy`

      Warning: If you have an older version of virtualenv you should update pip
      :code:`pip install --upgrade pip`.

      **Using Conda**

      Alternatively you can use conda for enviroment and package management. You
      can follow the installation instructions for windows users.

   .. tab-item:: Windows

      For windows we recommend using conda as package manager. You can download a
      light weight open source variant of conda: "miniforge3".

      1. Download latest `miniforge3 <https://github.com/conda-forge/miniforge>`__
         for Python 3.x (64 or 32 bit).
      2. Install miniforge3
      3. Open "miniforge prompt" to manage your virtual environments. You can
         create a new environment and acivate it by

         .. code-block:: console

            conda create -n tespy-env python=3.9
            activate tespy-env

      4. In the active prompt type: :code:`pip install tespy`
