import os
import sys

sys.path.insert(0, os.path.abspath("../src"))  # adjust path as needed
import soil_heat

# -- Project information -----------------------------------------------------
project = "Soil Heat"
copyright = "2025, Paul Inkenbrandt"
author = "Paul Inkenbrandt"
# The short X.Y version.
version = soil_heat.__version__
# The full version, including alpha/beta/rc tags.
release = soil_heat.__version__

extensions = [
    "sphinx.ext.autodoc",
    "numpydoc",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",
    "nbsphinx",
]

# Tell myst-parser to assign header anchors for h1-h3.
myst_heading_anchors = 3

suppress_warnings = ["myst.header"]

exclude_patterns = [
    "tests/*",
    "_build/*",
    "docs/_build/*",
    "Thumbs.db",
    ".DS_Store",
]  # Exclude the tests directory and _build directory

napoleon_numpy_docstring = True  # Set this to True for NumPy-style
napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
autosummary_generate = True  # Automatically generate .rst files for modules
autosummary_imported_members = True
bibtex_bibfiles = ["refs.bib"]  # Your BibTeX file(s)
bibtex_reference_style = "author_year"  # Use author-year style for citations
bibtex_default_style = "plain"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

source_suffix = [".rst", ".md"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    # Uncomment these if you use them in your codebase:
    #  "torch": ("https://pytorch.org/docs/stable", None),
    #  "datasets": ("https://huggingface.co/docs/datasets/master/en", None),
    #  "transformers": ("https://huggingface.co/docs/transformers/master/en", None),
}
