"""
====================================================================================================
Static Utilities (:mod:`SimPEG.electromagnetics.utils`)
====================================================================================================
.. currentmodule:: SimPEG.electromagnetics.static.utils


Electrode Utilities
===================
.. autosummary::
  :toctree: generated/

  electrode_separations
  pseudo_locations
  geometric_factor
  apparent_resistivity_from_voltage
  convert_survey_3d_to_2d_lines

Plotting Utilities
==================
.. autosummary::
  :toctree: generated/

  plot_pseudosection
  plot_3d_pseudosection


Survey Creation Utilities
=========================
.. autosummary::
  :toctree: generated/

  generate_survey_from_abmn_locations
  generate_dcip_survey
  generate_dcip_sources_line
  gen_3d_survey_from_2d_lines

Mesh Utilities
==============
.. autosummary::
  :toctree: generated/

  gettopoCC
  drapeTopotoLoc
  genTopography
  closestPointsGrid

"""
from .static_utils import (
    electrode_separations,
    pseudo_locations,
    geometric_factor,
    apparent_resistivity_from_voltage,
    plot_pseudosection,
    generate_dcip_survey,
    generate_dcip_sources_line,
    generate_survey_from_abmn_locations,
    convert_survey_3d_to_2d_lines,
    xy_2_lineID,
    r_unit,
    gettopoCC,
    drapeTopotoLoc,
    genTopography,
    closestPointsGrid,
    gen_3d_survey_from_2d_lines,
    plot_1d_layer_model,
)

# Import if user has plotly
try:
    from .static_utils import plot_3d_pseudosection
except ImportError:
    pass
