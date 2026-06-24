"""
cicoil_julia_hook.py — Alternative to the config-based integration.

For the recommended approach, set these in your run script:
    o.set_config('processes:julia_weathering', True)
    o.set_config('julia:project_path', '/path/to/CICOILPhysics.jl')

This file provides the monkey-patch approach for backward compatibility:
    from cicoil_julia_hook import patch_cicoil_julia
    patch_cicoil_julia(o, julia_project='path/to/CICOILPhysics.jl')

Both approaches are equivalent. The config method is preferred because it
integrates with CICOILv2's logging and timer infrastructure.
"""

import logging
import os

logger = logging.getLogger(__name__)


def patch_cicoil_julia(o, julia_project='CICOILPhysics.jl'):
    """Enable Julia weathering via CICOILv2 config (no monkey-patching).

    Parameters
    ----------
    o : OpenCiceseOil
        Simulation object (before calling run()).
    julia_project : str
        Path to the CICOILPhysics.jl directory containing Project.toml.
    """
    julia_project = os.path.abspath(julia_project)
    o.set_config('processes:julia_weathering', True)
    o.set_config('julia:project_path', julia_project)
    logger.info('Julia weathering enabled via config (project: %s)',
                julia_project)
