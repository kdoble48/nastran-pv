"""Sphinx configuration for the nastran-pv docs site (GitHub Pages)."""
import sys
from pathlib import Path

DOCS = Path(__file__).parent
sys.path.insert(0, str(DOCS / "_ext"))
sys.path.insert(0, str(DOCS.parent / "src"))  # autodoc imports nastran_pv from the checkout (api.md)

project = "nastran-pv"
author = "nastran-pv authors"
copyright = "2026, nastran-pv authors"
extensions = ["myst_parser", "sphinx_design", "sphinx.ext.autodoc", "npv_gallery"]
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3
exclude_patterns = ["_build", "design", "models", "README.md"]
# the Python API page: signatures as written in the source (registry=REGISTRY, not its repr)
autodoc_member_order = "bysource"
autodoc_preserve_defaults = True
html_theme = "pydata_sphinx_theme"
html_title = "nastran-pv"
html_static_path = ["_static"]
html_css_files = ["npv.css"]
html_theme_options = {
    "github_url": "https://github.com/kdoble48/nastran-pv",
    "navbar_align": "left",
    "header_links_before_dropdown": 6,
    "show_toc_level": 2,
}
