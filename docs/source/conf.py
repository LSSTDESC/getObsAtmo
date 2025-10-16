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

sys.path.insert(0, os.path.abspath("../.."))
import getObsAtmo  # noqa: E402


release = getObsAtmo._version.__version__  # version complète, ex: '0.2.1'
version = ".".join(release.split(".")[:2])  # version courte, ex: '0.2'


# sys.path.insert(0, os.path.abspath('../../getObsAtmo'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "getObsAtmo"
copyright = "2025, Sylvie Dagoret-Campagne"
author = "Sylvie Dagoret-Campagne"
# release = "0.2.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_tabs.tabs",
    "numpydoc",
    "nbsphinx",
    "sphinx.ext.graphviz",
    "sphinx.ext.inheritance_diagram",
    "myst_nb",
]


nbsphinx_execute = "always"
nbsphinx_allow_errors = True
source_suffix = [".rst"]


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
