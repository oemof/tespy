.. _developing_tespy_label:

Developing TESPy
================

Thus far, TESPy has been developed solely by Francesco Witte at the university of applied
sciences Flensburg. We hope that many people can make use of this project and hope it will be a
community driven project in the future, as users might demand special components or flexible
implementations of characteristics, custom equations, basically what ever you can think of.

We want to invite you to join the developing process and share your ideas. Your solutions may help
other users as well. Contributing to the development of TESPy may help other people and it improves
the quality of your code as it will be reviewd by other developers.

The easiest way of joining the developing process is by improving the documentation. If you spot
mistakes or if you think the documentation could be more precise or clear in some sections,
feel free to fix it and create a pull request on the github repository.

Another simple first step is to program your own custom subsystems and share them in the community. For
a good start just look into our tutorial `How do I create custom subsystems
<https://tespy.readthedocs.io/tutorial/custom_subsystem.html>`_.

A third, highly appreciated way for you to contribute is the provision of new and/or improved equations,
maybe even whole characteristics for single components, such as compressor maps, characteristics
for turbine efficiency or heat transfer coefficients etc..

You will find the most important information concerning the development process in the following
sections. If you have any further questions feel free to contact us, we are looking forward to hearing
from you!

.. contents::
    :depth: 1
    :local:
    :backlinks: top

Install the developer version
-----------------------------

It is recommenden to use `virtual environments <https://docs.python.org/3/tutorial/venv.html>`_.
Clone the github repository and install the repository in your virtual environment with pip:

.. code:: bash

  pip install -e /path/to/the/repository
   
If you already installed a previous developer version and want to update just use the following
within your environment.

.. code:: bash

  pip install --upgrade -e /path/to/the/repository
  
Contribute to the documentation
-------------------------------

If you come across typos or grammatical mistakes or want to improve comprehensibility of the
documentation, make your adjustments or suggestions and create a pull request for the dev branch.
We appreciate your contribution!

Improve or add new equations
----------------------------

The components equations represent the behaviour of each component. Are we missing any equations?
Or do you see possibilities to improve the formulation? Add them to your code and run some the
tests. If your tests are succesful create a pull request on the dev branch.

Share your subsystems
---------------------

You have encountered component groups that are frequently required? Make a subsystem to simplify
your own work and share your subsystem by adding it to the tespy.components.subsystems file.

For both, **the equations and the subsystems**, please provide the docstring and maybe inline comments,
if you want to make the code more comprehensible. In general, your code should fit to the
:ref:`style_guidlines_label`, a complete documentation is only required, if the code will be added
to TESPy permanently.

Add component characteristics
-----------------------------

The component characteristics represent large added value for your calculation. If you have detailed
information on components offdesign behaviour - even for specific cases - it will improve the results.
Every user can benefit from this knowlegde and thus we are very happy to discuss about the implementation
of new characteristics.

Collaboration with pull requests
--------------------------------

To collaborate use the pull request functionality of github as described here:
https://guides.github.com/activities/hello-world/

How to create a pull request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Fork the oemof repository to your own github account.
* Change, add or remove code.
* Commit your changes.
* Create a pull request and describe what you will do and why. Please use the pull request
  template we offer. It will be shown to you when you click on "New pull request".
* Wait for approval.

.. _coding_requirements_label:  

Generally the following steps are required when changing, adding or removing code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Read the :ref:`style_guidlines_label` and :ref:`naming_conventions_label` and follow them
* Add new tests according to what you have done
* Add/change the documentation (new feature, API changes ...)
* Add a whatsnew entry and your name to Contributors
* Check if all :ref:`tests_label` still work.

.. _tests_label:

Tests
-----

There are no tests available, yet. If you like to design a test, please contact us.

.. _style_guidlines_label:

Issue-Management
----------------

A good way for communication with the developer group are issues. If you
find a bug, want to contribute an enhancement or have a question on a specific problem
in development you want to discuss, please create an issue:

* describing your point accurately
* using the list of category tags
* addressing other developers

If you want to address other developers you can use @name-of-developer, or
use e.g. @tespy to address a team. `Here <https://github.com/orgs/oemof/teams>`_
you can find an overview over existing teams on different subjects and their members.

Look at the existing issues to get an idea on the usage of issues.

Style guidelines
----------------

We mostly follow standard guidelines instead of developing own rules. So if anything is
not defined in this section, search for a `PEP rule <https://www.python.org/dev/peps/>`_
and follow it.

Docstrings
^^^^^^^^^^

We decided to use the style of the numpydoc docstrings. See the following link for an
`example <https://github.com/numpy/numpy/blob/master/doc/example.py>`_.


Code commenting
^^^^^^^^^^^^^^^^

Code comments are block and inline comments in the source code. They can help to understand
the code and should be utilized "as much as necessary, as little as possible". When writing
comments follow the PEP 0008 style guide: https://www.python.org/dev/peps/pep-0008/#comments.


PEP8 (Python Style Guide)
^^^^^^^^^^^^^^^^^^^^^^^^^

* We adhere to `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ for any code
  produced in the framework.

* We use pylint to check your code. Pylint is integrated in many IDEs and 
  Editors. `Check here <http://docs.pylint.org/ide-integration>`_ or ask the 
  maintainer of your IDE or Editor

* Some IDEs have pep8 checkers, which are very helpful, especially for python 
  beginners.

Quoted strings
^^^^^^^^^^^^^^

As there is no recommendation in the PEP rules we use double quotes for strings read by humans
such as logging/error messages and single quotes for internal strings such as keys and column
names. However one can deviate from this rules if the string contains a double or single quote
to avoid escape characters. According to `PEP 257 <http://legacy.python.org/dev/peps/pep-0257/>`_
and numpydoc we use three double quotes for docstrings.

.. code-block:: python

    logging.info("We use double quotes for messages")
    
    my_dictionary.get('key_string')
    
    logging.warning('Use three " to quote docstrings!'  # exception to avoid escape characters

.. _naming_conventions_label:

Naming Conventions
------------------

* We use plural in the code for modules if there is possibly more than one child
  class (e.g. import transformers AND NOT transformer). If there are arrays in
  the code that contain multiple elements they have to be named in plural (e.g.
  `transformers = [T1, T2,...]`).

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
`Vincent Driessen <http://nvie.com/posts/a-successful-git-branching-model/>`_.

Differences are:

* instead of the name ``origin/develop`` we call the branch ``origin/dev``.
* feature branches are named like ``features/*``
* release branches are named like ``releases/*``

Commit message
^^^^^^^^^^^^^^

Use this nice little `commit tutorial <http://chris.beams.io/posts/git-commit/>`_ to 
learn how to write a nice commit message.
 

Documentation
----------------

The general implementation-independent documentation such as installation guide, flow charts,
and mathematical models is done via ReStructuredText (rst). The files can be found in the folder
*/oemof/doc*. For further information on restructured text see: http://docutils.sourceforge.net/rst.html.


