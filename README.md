# CICOILv2

**Component-resolved oil spill weathering model for the Gulf of Mexico**

CICOILv2 extends [OpenDrift](https://github.com/OpenDrift/opendrift) (v1.14.9) with physically-based, component-resolved oil weathering processes developed at CICESE (Centro de Investigacion Cientifica y de Educacion Superior de Ensenada). It tracks 17 pseudo-components through evaporation, emulsification, biodegradation, photooxidation, and subsea dissolution, with optional coupling to the [TAMOC](https://github.com/socolofs/tamoc) bent plume model for blowout scenarios.

## Features

- **17-component mass balance** — 8 aliphatic + 8 aromatic pseudo-groups + residue, resolved from ADIOS oil database distillation curves
- **Biodegradation** — Q10-scaled first-order decay with droplet-size dependence and dissolved-phase pathway
- **Photooxidation** — UV-driven aromatic degradation with Beer-Lambert depth attenuation; oxygenated photoproduct (OP) tracking with dissolved/degraded pools
- **Emulsification** — Fingas & Fieldhouse viscosity-stability model with configurable rate scaling
- **TAMOC coupling** — Bent plume model output seeded as Lagrangian elements
- **NEMO/CROCO readers** — Optimized NetCDF readers for regional ocean model output
- **97 element variables** tracked per particle

## Requirements

- Python 3.8+
- OpenDrift 1.14.9
- NumPy, SciPy, xarray, netCDF4, matplotlib
- ADIOS oil database (bundled with OpenDrift)
- Optional: TAMOC (for subsea blowout coupling)

## Installation

CICOILv2 operates as an in-place extension of OpenDrift. Install OpenDrift first, then place the CICOILv2 files in the OpenDrift package directory:

```bash
# 1. Create conda environment
conda create -n cicoil python=3.10
conda activate cicoil

# 2. Install OpenDrift 1.14.9
pip install opendrift==1.14.9 future   # 'future' provides the `past` module used by fluid_properties.py

# 3. Copy CICOILv2 files into OpenDrift's openoil module
SITE=$(python -c "import opendrift; print(opendrift.__path__[0])")
cp ciceseoil.py fluid_properties.py cicoil_estimations.py \
   tamoc_plume.py run_tamoc.py tamoc_chemical_properties.py \
   cicoil_common.py "$SITE/models/openoil/"
cp -r readers/* "$SITE/readers/"
mkdir -p "$SITE/models/openoil/data"
cp -r data/* "$SITE/models/openoil/data/"
cp export/io_stat_nc.py "$SITE/export/"
```

For HPC deployment (CHAMAN2), see `deploy_on_hpc_conda.sh`.

## Quick start

```python
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_constant import Reader as ConstReader
from datetime import datetime, timedelta

o = OpenCiceseOil(weathering_model='cicese')
o.set_oiltype('MAYA')

r = ConstReader({
    'x_sea_water_velocity': 0.10,
    'y_sea_water_velocity': 0.00,
    'x_wind': 5.0, 'y_wind': 0.0,
    'sea_water_temperature': 302.0,
    'sea_water_salinity': 36.0,
    'sea_floor_depth_below_sea_level': 50.0,
    'sea_surface_height': 0.0,
    'ocean_vertical_diffusivity': 0.01,
    'sea_surface_wave_significant_height': 1.0,
    'sea_surface_wave_stokes_drift_x_velocity': 0.02,
    'sea_surface_wave_stokes_drift_y_velocity': 0.0,
    'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.0,
    'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 5.0,
})
o.add_reader(r)

o.set_config('processes:biodegradation', True)
o.set_config('processes:photooxidation', True)

o.seed_elements(lon=-92.1, lat=19.7, number=200,
                time=datetime(2023, 7, 6, 12), m3_per_hour=100.)
o.run(duration=timedelta(hours=24), time_step=timedelta(hours=1))

budget = o.get_oil_budget()
```

## Directory structure

```
CICOILv2/
├── ciceseoil.py              # Main model class (OpenCiceseOil + CiceseOil element)
├── fluid_properties.py       # 17-component distillation framework
├── cicoil_estimations.py     # SARA fraction estimation
├── tamoc_plume.py            # TAMOC bent plume → OpenDrift seeding
├── run_tamoc.py              # TAMOC simulation wrapper
├── tamoc_chemical_properties.py
├── cicoil_common.py          # Shared helpers (load_ciceseoil, default_cicese_config)
├── plotting.py               # Publication plotting routines
├── download_data.py          # GLORYS12/ERA5 download + preprocessing
├── preprocess_era5_uv.py     # ERA5 UV radiation preprocessing
├── preprocess_real_forcing.py # Multi-source forcing assembly
├── readers/
│   ├── reader_nemo_optimized.py
│   ├── reader_NEMO_native_v3.py
│   ├── reader_nemo_combined.py
│   └── reader_nemo_modified.py
├── export/
│   └── io_stat_nc.py         # Gridded statistical NetCDF export
├── data/                     # Chemical property CSVs
├── test_cicoil_e2e.py        # 30-test end-to-end suite (T1–T8)
├── test_biodegradation.py
├── test_biodegradation_water.py
├── test_droplet_size_biodeg.py
├── test_emulsification_viscosity.py
├── test_op_tracking.py
├── test_photoox_era5uv.py
├── test_photoox_vertmix.py
├── deploy_on_hpc_conda.sh    # HPC deployment script
├── CICOIL_v2_MANUAL.md       # Full manual (physics + deployment)
├── CICOIL_v2_REFERENCE.md    # API reference
└── BIODEGRADATION.md         # Biodegradation module documentation
```

## Testing

```bash
conda run -n cicoil python test_cicoil_e2e.py
# Expected: PASSED: 30/30
```

The test suite covers:
- T1–T3: NOAA and CICESE weathering modes, mass conservation
- T4: TAMOC plume seeding
- T5: NEMO reader imports
- T6: Mass budget plotting
- T7: Biodegradation (component-level mass balance)
- T8: Photooxidation (OP tracking + regression bounds)

## Documentation

- **[CICOIL_v2_MANUAL.md](CICOIL_v2_MANUAL.md)** — Full manual: installation, deployment, weathering physics, parameter calibration
- **[CICOIL_v2_REFERENCE.md](CICOIL_v2_REFERENCE.md)** — API reference: element variables, configuration keys, method signatures

## Authors

- **Konstantinos Kotzakoulakis** (CICESE / SINTEF Ocean) — original CICOILv2 development
- **Knut-Frode Dagestad** (MET Norway) — OpenDrift framework
- **Alex Dominguez** (CICESE) — biodegradation, photooxidation, OP tracking, code corrections

## License

GNU General Public License v2.0 — see [LICENSE](LICENSE).

CICOILv2 is a derivative work of [OpenDrift](https://github.com/OpenDrift/opendrift), distributed under the same license.
