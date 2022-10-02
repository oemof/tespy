# -*- coding: utf-8 -*-

import os
import sys
import tespy


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('..'))

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

# landing page
# master_doc = 'contents'
# names, years, etc
project = 'TESPy'
year = '2022'
author = 'Francesco Witte'
copyright = '{0}, {1}'.format(year, author)

# The short X.Y version.
version = tespy.__version__.split(' ')[0]
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
pygments_style = 'sphinx'
pygments_style = 'trac'

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

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".


# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = '%s-%s' % (project, version)

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
}

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

# copybutton configuration
copybutton_prompt_text = r'>>> |\.\.\. '
copybutton_prompt_is_regexp = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'tespy_doc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'tespy.tex', u'tespy Documentation',
   u'Francesco Witte', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'tespy', u'tespy Documentation',
     [u'Francesco Witte'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'tespy', u'tespy Documentation',
   u'Author', 'tespy', 'Thermal Engineering Systems in Python',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False


# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = u'tespy'
epub_author = u'Francesco Witte'
epub_publisher = u'Francesco Witte'
epub_copyright = u'2018-2020, Francesco Witte'

# The basename for the epub file. It defaults to the project name.
#epub_basename = u'pahesmf'

# The HTML theme for the epub output. Since the default themes are not optimized
# for small screen space, using the same theme for HTML and epub output is
# usually not wise. This defaults to 'epub', a theme designed to save visual
# space.
#epub_theme = 'epub'

# The language of the text. It defaults to the language option
# or en if the language is not set.
#epub_language = ''

# The scheme of the identifier. Typical schemes are ISBN or URL.
#epub_scheme = ''

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#epub_identifier = ''

# A unique identification for the text.
#epub_uid = ''

# A tuple containing the cover image and cover page html template filenames.
#epub_cover = ()

# A sequence of (type, uri, title) tuples for the guide element of content.opf.
#epub_guide = ()

# HTML files that should be inserted before the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_pre_files = []

# HTML files shat should be inserted after the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_post_files = []

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# The depth of the table of contents in toc.ncx.
#epub_tocdepth = 3

# Allow duplicate toc entries.
#epub_tocdup = True

# Choose between 'default' and 'includehidden'.
#epub_tocscope = 'default'

# Fix unsupported image types using the PIL.
#epub_fix_images = False

# Scale large images.
#epub_max_image_width = 0

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#epub_show_urls = 'inline'

# If false, no index is generated.
#epub_use_index = True

# -- Options for linkcheck ----------------------------------------------

linkcheck_ignore = [
  r'https://doi.org/10.2172/95571',
]
