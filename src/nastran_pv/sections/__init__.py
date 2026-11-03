"""Beam cross sections: shape outlines, section properties and the 3D sweep.

* :mod:`.shapes` - outline rings of the PBARL/PBEAML library shapes (pure numpy).
* :mod:`.props` - exact polygon properties and the warping finite-element solve for J,
  the shear center and Cw; :func:`.props.library_section` / :func:`.props.equivalent_section`.
* :mod:`.beams` - the deck side: CBAR/CBEAM frames and offsets, property stations and the
  ``SECTION_*``/``SECRING_*``/``SECPT_*`` root tables.
* :mod:`.sweep` - sweeps the sections along the ``Lines`` block (no deck needed), rigid-section
  displacement, axes and triads.
"""
