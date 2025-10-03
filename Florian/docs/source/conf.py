# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Null Space Optimizer'
copyright = '2023, Florian Feppon'
author = 'Florian Feppon'
release = '1.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import sys  
import os
sys.path.append(os.path.abspath('lib'))
sys.path.insert(0,os.path.abspath('../..'))
import nullspace_optimizer 
    
#extensions = []
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary', 
              'sphinxawesome_theme',    
              'sphinx_contrib_dir',  
              'sphinx_subfigure',
              'sphinxcontrib.katex']

templates_path = ['_templates']
exclude_patterns = []

katex_prerender = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_extra_path = ["google7bdbeab653ef60e9.html"]
html_theme = "sphinxawesome_theme"
    
add_module_names = False
autodoc_member_order = 'bysource'
