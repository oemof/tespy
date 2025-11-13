# -*- coding: utf-8 -*-

import os
import sys
import tespy
from sphinx.ext import apidoc


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

THIS_PATH = os.path.dirname(__file__)
API_PATH = os.path.join(THIS_PATH, "api")
TESPY_PATH = os.path.join(THIS_PATH, "..", "src", "tespy")
# -- General configuration ------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
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


def modify_rst_contents():
    pass


def _check_rst_file_empty():
    return False


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


apidoc.main(["-f", "-M", "-d", "1", "-o", API_PATH, TESPY_PATH])

source_files = [file for file in os.listdir(API_PATH) if file.endswith(".rst")]
for file in source_files:
    _remove_subpackages_list_and_submodules_header(os.path.join(API_PATH, file))


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
exclude_patterns = ['_build']

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
