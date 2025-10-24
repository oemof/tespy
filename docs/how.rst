.. _tespy_development_how_label:

How can I contribute?
=====================

You will find the most important information concerning the development process
in the following sections. If you have any further questions feel free to
contact us, we are looking forward to hearing from you!

Install the developer version
-----------------------------

It is recommenden to use
`virtual environments <https://docs.python.org/3/tutorial/venv.html>`_ for
the development process. Fork the repository and clone your forked tespy GitHub
repository and install development requirements with pip.

.. code:: bash

    git clone https://github.com/YOUR_GITHUB_USERNAME/tespy.git
    cd tespy
    pip install -e .[dev]

In order to stay in sync with the oemof/tespy base repository, add the link to
the oemof/tespy repository as remote to your local copy of tespy. We will call
the link "upstream". Fetch branches available and after that, you can pull
changes from a specific branch of the oemof/tespy repository (branch "dev" in
the example below).

.. code:: bash

    git remote add upstream https://github.com/oemof/tespy.git
    git fetch upstream
    git pull upstream dev --rebase

Use the :code:`--rebase` command to avoid merge commits for every upstream pull.
If you want to make changes to tespy, checkout a new branch from your local dev
branch. Make your changes, commit them and create a PR on the oemof/tespy dev
branch.

Collaboration with pull requests
--------------------------------

To collaborate use the pull request functionality of GitHub as described here:
https://guides.github.com/activities/hello-world/

How to create a pull request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Fork the oemof repository to your own GitHub account.
* Change, add or remove code.
* Commit your changes.
* Create a pull request and describe what you will do and why. Please use the
  pull request template we offer. It will be shown to you when you click on
  "New pull request".
* Wait for approval.

.. _coding_requirements_label:

Generally the following steps are required when changing, adding or removing code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Read the :ref:`style_guidlines_label` and :ref:`naming_conventions_label` and
  follow them
* Add new tests according to what you have done
* Add/change the documentation (new feature, API changes ...)
* Add a What's New entry and your name to the list of contributors
* Check if all :ref:`tests_label` still work.

.. _tests_label:

Tests
-----

The tests in TESPy are split up in two different parts:

* doc-tests (also used as examples for classes and methods/functions)
* software tests (defined in the tests folder).

The tests contain code examples that expect a certain outcome. If the outcome
is as expected a test will pass, if the outcome is different, the test will
fail. You can run the tests locally by navigating into your local GitHub clone.
The command :code:`check` tests PEP guidelines, the command :code:`docs`
tests building the documentation, and the command :code:`py3X` runs the
software tests in the selected Python version.

.. code:: bash

    python -m tox -e docs
    python -m tox -e check
    python -m tox -e py36
    python -m tox -e py37
    python -m tox -e py38

If you want to have a look at the documentation build on your local machine use
the following command from the local tespy clone:

.. code:: bash

    python -m sphinx docs/ path/to/html_output

Additionally, all tests will run automatically when you push changes to a
branch that has a pull request opened.

If you have further questions regarding the tests, we are looking forward to
your inquiry.

.. _style_guidlines_label:

Issue-Management
----------------

A good way for communication with the developer group are issues. If you
find a bug, want to contribute an enhancement or have a question on a specific
problem in development you want to discuss, please create an issue:

* describing your point accurately
* using the list of category tags
* addressing other developers

If you want to address other developers you can use @name-of-developer, or
use e.g. @tespy to address a team.
`Here <https://github.com/orgs/oemof/teams>`__ you can find an overview over
existing teams on different subjects and their members.

Look at the existing issues to get an idea on the usage of issues.

Style guidelines
----------------

We mostly follow standard guidelines instead of developing own rules. So if
anything is not defined in this section, search for a
`PEP rule <https://www.python.org/dev/peps/>`_ and follow it.

Docstrings
^^^^^^^^^^

We decided to use the style of the numpydoc docstrings. See the following
link for more information
`numpy docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.


Code commenting
^^^^^^^^^^^^^^^^

Code comments are block and inline comments in the source code. They can help
to understand he code and should be utilized "as much as necessary, as little
as possible". When writing comments follow the
`PEP 0008 style guide <https://www.python.org/dev/peps/pep-0008/#comments>`_.


PEP8 (Python Style Guide)
^^^^^^^^^^^^^^^^^^^^^^^^^

* We adhere to `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ for any code
  produced in the framework.

* We use pylint to check your code. Pylint is integrated in many IDEs and
  Editors. `Check here <https://pylint.pycqa.org/en/latest/>`_ or ask the
  maintainer of your IDE or Editor

* Some IDEs have pep8 checkers, which are very helpful, especially for python
  beginners.

.. _naming_conventions_label:

Naming Conventions
------------------

* We use plural in the code for modules if there is possibly more than one
  child class (e.g. :code:`import heat_exchangers` AND NOT
  :code:`import heat_exchanger`). If there are arrays in the code that contain
  multiple elements they have to be named in plural.

* Please, follow the naming conventions of
  `pylint <http://pylint-messages.wikidot.com/messages:c0103>`_

* Use talking names

  * Variables/Objects: Name it after the data they describe
    (power\_line, wind\_speed)
  * Functions/Method: Name it after what they do: **use verbs**
    (get\_wind\_speed, set\_parameter)


Using git
---------

Branching model
^^^^^^^^^^^^^^^

So far we adhere mostly to the git branching model by
`Vincent Driessen <https://nvie.com/posts/a-successful-git-branching-model/>`_.

Differences are:

* instead of the name ``origin/develop`` we call the branch ``origin/dev``.
* feature branches are named like ``features/*``
* release branches are named like ``releases/*``

Commit message
^^^^^^^^^^^^^^

Use this nice little `commit tutorial <https://commit.style/>`_ to
learn how to write a nice commit message.


Documentation
-------------

The general implementation-independent documentation such as installation
guide, flow charts, and mathematical models is done via ReStructuredText (rst).
The files can be found in the folder *docs*. For further information on
restructured text see: https://docutils.sourceforge.io/rst.html.


Examples
--------

.. toctree::
    :maxdepth: 1
    :hidden:

    dev_examples/new_comp.rst

.. card::  How to contribute a new component
    :link: tespy_development_new_comp_label
    :link-type: ref
    
    .. figure:: /_static/images/tutorials/compressor_with_cooling/drawing_pcwc.svg
        :align: center
        :alt: Blackboard drawing of compressor with cooling
        :figclass: only-light

    .. figure:: /_static/images/tutorials/compressor_with_cooling/drawing_pcwc_darkmode.svg
        :align: center
        :alt: Blackboard drawing of compressor with cooling
        :figclass: only-dark
