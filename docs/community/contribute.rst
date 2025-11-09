.. _tespy_development_what_label:

~~~~~~~~~~~~~~~~~~~
Contribute to TESPy
~~~~~~~~~~~~~~~~~~~

What can I contribute?
======================
TESPy has been developed mainly by Francesco Witte at the University of Applied
Sciences Flensburg - and since 2021 with a lot of free time. The goal is, that
many people can make use of this project and that it will be a community driven
project in the future, as users might demand special components or flexible
implementations of characteristics, custom equations, basically what ever you
can think of.

Therefore, I would like to invite you to contribute in this process, share your
ideas and experience and maybe start developing the software. Your solutions
may help other users as well. Contributing to the development of TESPy is easy
and will help the development team and all other users of the software. If you
start to develop features it may also improve the quality of your code as it
will be reviewed by other developers.

There are a variety of different ways you can contribute to the development,
amongst others, these may be:

Contribute to the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The easiest way of joining the developing process is by improving the
documentation. If you spot mistakes or think, the documentation could be more
precise or clear in some sections, feel free to fix it and create a pull
request on the GitHub repository.

If you come across typos or grammatical mistakes or want to improve
comprehensibility of the documentation, make your adjustments or suggestions
and create a pull request for the dev branch. We appreciate your contribution!

Share your projects
^^^^^^^^^^^^^^^^^^^
You have used the software in your research paper or project, maybe even in a
real world application? We would love to feature your project on our
:ref:`Example Applications <tespy_examples_label>` page. Please reach out to
us by opening a new issue on our GitHub page.

Add new component equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The components equations represent the behavior of each component. Do you miss
equations? Open a discussion on the GitHub discussions page or add them to your
fork of TESPy and create a pull request on the dev branch.

Sponsor the development
^^^^^^^^^^^^^^^^^^^^^^^
As mentioned above, I maintain and develop this project in my free time. If you
would like to support me by sponsoring a coffee for late night sessions or need
direct support for your projects using TESPy, you sponsor me on GitHub
https://github.com/fwitte or reach out to me via e-mail.

Add component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Another highly appreciated way for you to contribute is the provision of new
and/or improved characteristics for single components, such as compressor
efficiency maps, characteristics for turbine efficiency or heat transfer
coefficients etc..

The component characteristics represent large added value for your calculation.
If you have detailed information on components offdesign behavior - even for
specific cases - it will improve the results. Every user can benefit from this
knowledge, and thus we are very happy to discuss the implementation of new
characteristics.

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
    python -m tox -e py310
    python -m tox -e py311
    python -m tox -e py312

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
----------------

The general implementation-independent documentation such as installation
guide, flow charts, and mathematical models is done via ReStructuredText (rst).
The files can be found in the folder *docs*. For further information on
restructured text see: https://docutils.sourceforge.io/rst.html.
