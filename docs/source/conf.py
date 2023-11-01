# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import pathlib
sys.path.insert(0, os.path.abspath('../../getObsAtmo'))
#sys.path.insert(0, os.path.abspath('../../'))
#sys.path.insert(0,'/Users/dagoret/MacOSX/GitHub/LSST/getObsAtmo/getObsAtmo/')
#sys.path.insert(0, os.path.abspath("."))
#sys.path.insert(0, os.path.abspath("../"))
#sys.path.insert(0, os.path.abspath("../../getObsAtmo"))
#sys.path.insert(1, os.path.dirname(os.path.abspath("../")) + os.sep + "feature_engine")

#sys.path.append('/Users/dagoret/MacOSX/GitHub/LSST/getObsAtmo/docs/notebooks')
#sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())



# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'getObsAtmo'
copyright = '2023, Sylvie Dagoret-Campagne'
author = 'Sylvie Dagoret-Campagne'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ 'sphinx.ext.autodoc',
               'sphinx.ext.duration',
               'sphinx.ext.doctest',
               'sphinx.ext.intersphinx',
               'sphinx.ext.autosummary',
               'sphinx.ext.mathjax', 
               'sphinx.ext.napoleon', 
               'sphinx.ext.viewcode',
               'sphinx_tabs.tabs',
               'numpydoc',
               'nbsphinx',
               'sphinx.ext.graphviz',
               'sphinx.ext.inheritance_diagram',
               ]



nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
source_suffix = ['.rst']



templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']



