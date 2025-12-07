# -*- coding: utf-8 -*-

import os

import pandas as pd
from sphinx.ext import apidoc
from tabulate import tabulate

import tespy
from tespy.components import Source, Sink, PowerSource, PowerSink
from tespy.components.component import component_registry
from tespy.connections.connection import connection_registry
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import ReferencedFluidProperties as dc_ref
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

DOCS_ROOT_PATH = os.path.join(os.path.dirname(__file__))
API_PATH = os.path.join(DOCS_ROOT_PATH, "api")
TESPY_PATH = os.path.join(DOCS_ROOT_PATH, "..",  "src", "tespy")


def _create_component_instances() -> dict:
    instances = {}
    for cls_name, cls_ref in component_registry.items.items():
        instances[cls_name] = cls_ref(cls_name)

    return instances


def _create_connection_instances() -> dict:
    instances = {}
    for cls_name, cls_ref in connection_registry.items.items():
        if cls_name == "Connection":
            instance = cls_ref(Source("source"), "out1", Sink("sink"), "in1")
        elif cls_name == "PowerConnection":
            instance = cls_ref(PowerSource("source"), "power", PowerSink("sink"), "power")
        else:
            # do not add to the dict of not any of both classes
            continue

        instances[cls_name] = instance

    return instances


def _get_eq_reference(datacontainer):
    if not hasattr(datacontainer, "func") or datacontainer.func is None and datacontainer.structure_matrix is None:
        return

    elif datacontainer.func is None:
        basename = datacontainer.structure_matrix.__qualname__
        eq_reference = f"{datacontainer.structure_matrix.__module__}.{basename}"

    else:
        basename = datacontainer.func.__qualname__
        eq_reference = f"{datacontainer.func.__module__}.{basename}"

    return ":py:meth:`" + basename.split(".")[-1] + " <" + eq_reference + ">`"


def collect_component_parameters(instance):
    grouped_parameters = {}
    parameters_with_equation = {}
    characteristic_lines_and_maps = {}
    constraints = {}

    for key, value in instance.get_mandatory_constraints().items():
        eq_reference = _get_eq_reference(value)
        constraints[key] = {
            "eq_reference": eq_reference,
            "description": value.description
        }

    for key, value in instance.parameters.items():
        eq_reference = _get_eq_reference(value)

        if isinstance(value, dc_cp):
            if value._potential_var:
                key = key + " [1]_"
            parameters_with_equation[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity": value.quantity
            }
        elif isinstance(value, dc_simple):
            parameters_with_equation[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity":  None
            }
        elif isinstance(value, dc_cc) or isinstance(value, dc_cm):
            characteristic_lines_and_maps[key] = {
                "eq_reference": eq_reference,
                "description": value.description
            }
        # gcp and gcc are equivalent
        elif isinstance(value, dc_gcp):
            grouped_parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "elements": ", ".join([f":code:`{element}`" for element in value.elements])
            }

    return constraints, parameters_with_equation, grouped_parameters, characteristic_lines_and_maps


def collect_connection_parameters(instance):
    parameters = {}

    for key, value in instance.property_data.items():
        eq_reference = _get_eq_reference(value)

        if isinstance(value, dc_prop):
            parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity": value.quantity
            }
        elif isinstance(value, dc_flu):
            parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity": value.quantity
            }
        elif isinstance(value, dc_ref):
            parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity": value.quantity
            }
        elif isinstance(value, dc_simple):
            parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity":  None
            }

    return parameters


def collect_component_constraints():
    return


def _indent_block(s: str, spaces: int) -> str:
    pad = " " * spaces
    return s.replace("\n", "\n" + pad)


def create_tabular_component_views():
    instances = _create_component_instances()
    path = os.path.join(DOCS_ROOT_PATH, "building_blocks", "_components_overview.rst")

    modules = {}

    for cls_name, instance in instances.items():
        # this is required to collect the mandatory constraints
        # it will keep the power balance equation between the component and the
        # (optionally) attached PowerConnection out of the list of equations
        # TODO: Find a better solution
        instance.power_outl = []
        instance.power_inl = []

        cls_ref = instance.__class__
        cls_module = cls_ref.__module__
        parent_module = ".".join(cls_ref.__module__.split(".")[:-1])
        if parent_module not in modules:
            modules[parent_module] = {}

        constraints, parameters, parameter_groups, characteristic_lines = collect_component_parameters(instance)

        modules[parent_module][cls_name] = _indent_block("\n", 4)
        modules[parent_module][cls_name] += _indent_block(f".. dropdown:: {cls_name}" + "\n" * 2, 8)
        modules[parent_module][cls_name] += _indent_block(
            f"Class documentation and example: :py:class:`{cls_name} <{cls_module}.{cls_name}>`", 8
        )

        outputs = {
            "Table of constraints": constraints,
            "Table of parameters": parameters,
            "Table of parameter groups": parameter_groups,
            "Table of characteristic lines and maps": characteristic_lines,
        }
        headers = {
            "description": "Description",
            "quantity": "Quantity",
            "elements": "Required parameters",
            "eq_reference": "Method"
        }

        for key, value in outputs.items():
            df = pd.DataFrame.from_dict(value).T
            if not df.empty:
                df = df.fillna(":code:`None`")

                modules[parent_module][cls_name] += _indent_block("\n" * 2, 8)
                modules[parent_module][cls_name] += _indent_block(
                    key + "\n" * 2, 8
                )

                modules[parent_module][cls_name] += _indent_block(
                    tabulate(
                        df[[c for c in headers if c in df.columns]],
                        headers=["Parameter"] + [v for c, v in headers.items() if c in df.columns],
                        tablefmt="rst"
                    ), 8
                )

        modules[parent_module][cls_name] += ("\n" * 2)

    with open(path, "w", encoding="utf-8") as f:
        f.write(".. tab-set::\n\n")
        for module, data in modules.items():
            f.write(f"    .. tab-item:: {module.split('.')[-1]}" + "\n" * 2)
            f.write(f"        .. container:: accordion-group" + "\n")

            for cls_info in data.values():
                f.write(_indent_block(cls_info, 8))
                f.write("\n" * 2)

        f.write(
            ".. [1] This parameter can be made a variable, "
            ":ref:`get more info here <component_variables_label>`."
        )


def create_tabular_connection_views():
    instances = _create_connection_instances()
    path = os.path.join(DOCS_ROOT_PATH, "building_blocks", "_connections_overview.rst")

    classes = {}

    for cls_name, instance in instances.items():

        cls_ref = instance.__class__
        cls_module = cls_ref.__module__

        parameters = collect_connection_parameters(instance)

        classes[cls_name] = "\n"
        classes[cls_name] += (f"Class documentation and example: :py:class:`{cls_name} <{cls_module}.{cls_name}>`")

        outputs = {
            "Table of parameters": parameters
        }
        headers = {
            "description": "Description",
            "quantity": "Quantity",
            "eq_reference": "Method"
        }

        for key, value in outputs.items():
            df = pd.DataFrame.from_dict(value).T
            if not df.empty:
                df = df.fillna(":code:`None`")

                classes[cls_name] += ("\n" * 2)
                classes[cls_name] += (key + "\n" * 2)
                classes[cls_name] += (
                    tabulate(
                        df[[c for c in headers if c in df.columns]],
                        headers=["Parameter"] + [v for c, v in headers.items() if c in df.columns],
                        tablefmt="rst"
                    )
                )

        classes[cls_name] += ("\n" * 2)

    with open(path, "w", encoding="utf-8") as f:
        f.write(".. tab-set::\n\n")
        for cls_name, cls_info in classes.items():
            f.write(f"    .. tab-item:: {cls_name.split('.')[-1]}\n")
            f.write(_indent_block(cls_info, 8))
            f.write("\n" * 2)


def _remove_subpackages_list_and_submodules_header(file):

    with open(file, "r") as f:
        lines = f.readlines()

    new_lines = []
    in_subpackages_section = False
    submodules_line = -2  # start somewhere guaranteed outside of range
    has_submodules = False
    for i, line in enumerate(lines):

        if line.strip() == "Subpackages":
            in_subpackages_section = True
            continue

        elif line.strip() == "Submodules":
            in_subpackages_section = False
            submodules_line = i
            has_submodules = True
            continue

        if not in_subpackages_section and i != submodules_line + 1:
            new_lines.append(line)

    if not has_submodules and "data.rst" not in file:
        os.remove(file)
    else:
        with open(file, "w") as f:
            f.write("".join(new_lines))


def update_api_docs():
    apidoc.main(["-f", "-M", "-d", "1", "-o", API_PATH, TESPY_PATH])

    source_files = [file for file in os.listdir(API_PATH) if file.endswith(".rst")]
    for file in source_files:
        _remove_subpackages_list_and_submodules_header(os.path.join(API_PATH, file))


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# -- General configuration ------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'myst_nb',
    'notfound.extension',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinx_design',
    'sphinxcontrib.bibtex',
]


def generate_all_sources(app):
    update_api_docs()
    create_tabular_component_views()
    create_tabular_connection_views()

def setup(app):
    app.connect("builder-inited", generate_all_sources)


# landing page
# master_doc = 'contents'
# names, years, etc
project = 'TESPy'
year = '2024'
author = 'Francesco Witte'
copyright = '{0}, {1}'.format(year, author)

# The short X.Y version.
version = tespy.__version__.split(' - ')[0]
# The full version, including alpha/beta/rc tags.
release = tespy.__version__

# The suffix of source filenames.
source_suffix = '.rst'
# folder for templates
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'jupyter_execute']

# The name of the Pygments (syntax highlighting) style to use.
# pygments_style = "some"
# pygments_dark_style = "someother"

# show all class members
# numpydoc_show_class_members = False

# place for bibtex references
bibtex_bibfiles = ['references.bib']

# links to github
github_repo_url = "https://github.com/oemof/tespy/"
extlinks = {
    "issue": (f"{github_repo_url}/issues/%s", "#%s"),  # noqa: WPS323
    "pr": (f"{github_repo_url}/pull/%s", "PR #%s"),  # noqa: WPS323
    "commit": (f"{github_repo_url}/commit/%s", "%s"),  # noqa: WPS323
}

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'furo'

html_title = f"{project} v{version}"

# Some more stuff
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
# html_additional_pages = {
#     "index": "index.html"
# }

html_sidebars = {
    '**': [
        'sidebar/brand.html',
        'sidebar/search.html',
        'sidebar/scroll-start.html',
        'sidebar/navigation.html',
        'sidebar/ethical-ads.html',
        'sidebar/scroll-end.html',
        'sidebar/variant-selector.html',
    ],
}

html_theme_options = {
    "light_logo": "./images/logo_tespy_mid.svg",
    "dark_logo": "./images/logo_tespy_mid_darkmode.svg",
    "announcement": """
    <div oemof-announcement=\"https://raw.githubusercontent.com/oemof/tespy/announcements/announcement.html\"></div>
    """,
    "light_css_variables": {
        "color-link": "#1f567d",
        "color-sidebar-link-text--top-level": "#1f567d",
    },
    "dark_css_variables": {
        "color-link": "#4a7ca3",
        "color-sidebar-link-text--top-level": "#4a7ca3",
    }
}

html_js_files = [
    'js/custom.js',
]


html_favicon = './_static/images/logo_tespy_small.svg'

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

# copybutton configuration
copybutton_prompt_text = r'>>> |\.\.\. '
copybutton_prompt_is_regexp = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'tespy_doc'

# -- Options for linkcheck ----------------------------------------------

linkcheck_ignore = [    # DOIs always redirect, we believe they will always work.
  r"https://doi.org/*",
  r'https://github\.com/oemof/tespy/.*',
  r'https://coolprop\.org/*'
]

# Notebook execution
nb_execution_allow_errors = True
nb_execution_in_temp = True
