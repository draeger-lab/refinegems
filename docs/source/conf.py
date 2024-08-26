# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path
sys.path.insert(0, os.path.abspath(Path('..','..','src')))

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True

# -- Project information -----------------------------------------------------

project = 'refineGEMs'
project_copyright = '2024, Famke Bäuerle, Gwendolyn O. Döbel and Carolin Brune'
author = 'Famke Bäuerle, Gwendolyn O. Döbel and Carolin Brune'

# The full version, including alpha/beta/rc tags
release = '2.0.0-alpha.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'nbsphinx',
    'sphinx_rtd_theme',
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinxcontrib.bibtex',
    'sphinx.ext.viewcode'
]

# For copy buttons in code blocks
copybutton_selector =  "div.copyable pre"

# For citations
bibtex_bibfiles = ['library.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme' #'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Change colours in theme for navigation
html_css_files = ['custom_theme.css']

# Adds logo to documentation page
html_logo = str(Path('images','refineGEMs_logo.png'))
html_theme_options = {
    'logo_only': True,
    'display_version': False,
	'flyout_display': 'attached'
}

#Adds logo as favicon to tab
html_favicon = str(Path('images','refineGEMs_logo.png'))

# Changes code highlighting
pygments_style = 'blinds-light'

# Make figures numbered
numfig = True

# Explicitly assign the master document
master_doc = 'index'

# -- Autodoc -----------------------------------------------------------------

autodoc_preserve_defaults = True
autodoc_default_options = {
    'exclude-members': '''
    __init__, __abstractmethods__, ABC, __module__, _abc_impl, __annotations__, 
    __dict__, __weakref__
    '''
}

# we need those to display the code comments otherwise the functions cannot be imported
# autodoc_mock_imports = ["psycopg2", 
#                         "gffutils",
#                         "cplex.exceptions",
#                         "cplex",
#                         "cobra",
#                         "pandas",
#                         "libsbml",
#                         "numpy",
#                         "bioservices",
#                         "bioregistry",
#                         "bs4",
#                         "memote",
#                         "tqdm",
#                         "psycopg2",
#                         "Bio",
#                         "sqlalchemy",
#                         "ratelimit",
#                         "libchebipy",
#                         "ols_client",
#                         "charges",
#                         "click",
#                         "yaml",
#                         "sortedcontainers",
#                         "colorama",
#                         "matplotlib",
#                         "seaborn",
#                         "venn"
#                         ]