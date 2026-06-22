"""
cicoil_common.py — shared helpers for CICOILv2 test and run scripts.

Eliminates the repeated importlib boilerplate and common set_config calls.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_ciceseoil(module_path=None):
    """Load OpenCiceseOil from a local ciceseoil.py without installing.

    Parameters
    ----------
    module_path : Path or str, optional
        Path to ciceseoil.py.  Defaults to the copy next to this file.

    Returns
    -------
    OpenCiceseOil : class
        The loaded model class, ready to instantiate.
    """
    if module_path is None:
        module_path = HERE / 'ciceseoil.py'
    module_path = Path(module_path)

    import opendrift.models.openoil  # noqa: E402
    spec = importlib.util.spec_from_file_location(
        'opendrift.models.openoil.ciceseoil_local', module_path)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = 'opendrift.models.openoil'
    sys.modules['opendrift.models.openoil.ciceseoil_local'] = mod
    spec.loader.exec_module(mod)

    print(f'Loaded OpenCiceseOil from: {module_path}')
    return mod.OpenCiceseOil


def default_cicese_config(o, **overrides):
    """Apply standard CICESE-mode config to an OpenCiceseOil instance.

    Default: evaporation+spreading+emulsification ON, everything else OFF,
    no vertical mixing.  Pass keyword overrides to change individual settings.

    Parameters
    ----------
    o : OpenCiceseOil instance
    **overrides : bool
        Keys are config suffixes (e.g. biodegradation=True).
        Supported keys: vertical_mixing, evaporation, spreading,
        emulsification, dispersion, biodegradation, photooxidation,
        handle_released_gas, subsea_dissolution.
    """
    defaults = {
        'vertical_mixing': False,
        'evaporation': True,
        'spreading': True,
        'emulsification': True,
        'dispersion': False,
        'biodegradation': False,
        'photooxidation': False,
        'handle_released_gas': False,
        'subsea_dissolution': False,
    }
    defaults.update(overrides)

    o.set_config('drift:vertical_mixing', defaults.pop('vertical_mixing'))
    for key, val in defaults.items():
        o.set_config(f'processes:{key}', val)
