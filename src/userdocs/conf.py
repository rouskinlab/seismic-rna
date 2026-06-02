# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Project information -----------------------------------------------------

project = "seismic-rna"
project_copyright = "2022-2026, the Rouskin Lab"
author = "Matthew Allan, Yves Martin, Scott Grote, Alberic de Lajarte, and Justin Aruda"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.

source_suffix = [".rst"]

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",
    "sphinx_click",
]

# Show todo directives in rendered output so reviewers see what still
# needs to be filled in.
todo_include_todos = True

# Templates
templates_path = ["_templates"]

# Static files (CSS, JS, images)
html_static_path = ["_static", "../../logo"]

# Generate automatic summaries
autosummary_generate = True

# Exclude from autosummary
exclude_patterns = ["_build", "_templates"]

# Fontpath for blockdiag (truetype font)
blockdiag_fontpath = "/usr/share/fonts/truetype/ipafont/ipagp.ttf"

# Provide a GitHub API token:
# Pass the SPHINX_GITHUB_CHANGELOG_TOKEN environment variable to your build
# OR
sphinx_github_changelog_token = "..."

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

# -- Options for HTML output
html_theme = "renku"
html_logo = "../../logo/logo-200.png"
html_favicon = "../../logo/favicon.ico"
html_css_files = ["seismic.css"]
html_js_files = ["nav_labels.js"]

html_theme_options = {
    "style_nav_header_background": "#a00821",
    "collapse_navigation": False,
}

# -- Options for EPUB output
epub_show_urls = "footnote"

autosectionlabel_prefix_document = True

# Options for build
# fail_on_warning = True
# autodoc_mock_imports = MOCK_MODULES
nitpick_ignore = [("py:class", "type")]
