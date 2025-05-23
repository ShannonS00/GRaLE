# Configuration file for the Sphinx documentation builder.

import os
import sys

#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

sys.path.insert(0, os.path.abspath('../../grale')) # path to file

# Add root dir (or examples dir) so Sphinx/nbsphinx can access notebooks
#sys.path.insert(0, os.path.abspath('../../Examples'))

project = 'GRaLE'
copyright = '2025, Shannon Schroeder'
author = 'Shannon Schroeder'
release = '1.0.0'
version = '1.0'  
language = 'en'


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'nbsphinx',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


html_theme = 'sphinx_rtd_theme'

# Include the static path for custom assets
html_static_path = ["_static"]

# Set the logo (appears top-left)
html_logo = "_static/logo.png"

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "ShannonS00",     # GitHub username
    "github_repo": "GRaLE",    # Repository name
    "github_version": "main",           # Branch
    "conf_py_path": "/docs/source/",    # Path in repo to your conf.py
}