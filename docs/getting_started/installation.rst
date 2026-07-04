.. _installation_and_setup_label:

######################
Installation and setup
######################

Following you find guidelines for the installation process for linux and
windows. TESPy is a Python package, thus it requires you to have Python 3
installed.

.. tab-set::

   .. tab-item:: uv (recommended)

      `uv <https://docs.astral.sh/uv/>`_ is a fast Python package and project
      manager. Install it by following the
      `official instructions <https://docs.astral.sh/uv/getting-started/installation/>`_,
      then initialise a new project and add TESPy:

      .. code-block:: console

         uv init my-project
         cd my-project
         uv add tespy

      To include the optional optimisation dependencies:

      .. code-block:: console

         uv add tespy[opt]

   .. tab-item:: pip

      We recommend installing TESPy within a virtual Python environment. Create
      and activate one, then install TESPy via pip:

      .. code-block:: console

         python -m venv tespy-env
         source tespy-env/bin/activate  # Linux/macOS
         # tespy-env\Scripts\activate   # Windows
         pip install tespy

      To include the optional optimisation dependencies:

      .. code-block:: console

         pip install tespy[opt]

   .. tab-item:: conda

      TESPy is available on conda-forge. You can install it into a new or
      existing conda environment:

      .. code-block:: console

         conda install -c conda-forge tespy

   .. tab-item:: Developer Version

      If you would like to get access to not yet released features or features
      under development you can install the developer version. Fork and clone the
      `TESPy repository <https://github.com/oemof/tespy>`_, then install all
      development dependencies with uv:

      .. code-block:: console

         cd tespy
         uv sync --extra dev

      Alternatively, follow the instructions on
      :ref:`this page <development_how_label>`.
