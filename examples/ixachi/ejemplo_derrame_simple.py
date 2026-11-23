#!/usr/bin/env python3
"""Ejemplo minimo de CICOILv2: derrame de MAYA en la Bahia de Campeche con forzantes constantes, 24 h.

No necesita datos externos (usa un lector constante), sirve para comprobar que el entorno funciona
y para ver las opciones de configuracion nuevas de v2.  NO es un caso fisico realista.
Uso:  python ejemplo_derrame_simple.py --out salida_ejemplo.nc [--n 200] [--dispersion-mode legacy|transport|budget]
"""
import argparse
from datetime import datetime, timedelta
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_constant import Reader as ConstReader

ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True)
ap.add_argument('--n', type=int, default=200)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--dispersion-mode', choices=['legacy', 'transport', 'budget'], default='legacy')
a = ap.parse_args()

# `seed=` va al constructor: OpenDrift controla su propio generador aleatorio
o = OpenCiceseOil(weathering_model='cicese', loglevel=30, seed=a.seed)
o.set_oiltype('MAYA')

o.add_reader(ConstReader({
    'x_sea_water_velocity': 0.10, 'y_sea_water_velocity': 0.0,
    'x_wind': 5.0, 'y_wind': 0.0,
    'sea_water_temperature': 302.0, 'sea_water_salinity': 36.0,
    'sea_floor_depth_below_sea_level': 50.0, 'sea_surface_height': 0.0,
    'ocean_vertical_diffusivity': 0.01,
    'sea_surface_wave_significant_height': 1.0,
    'sea_surface_wave_stokes_drift_x_velocity': 0.02,
    'sea_surface_wave_stokes_drift_y_velocity': 0.0,
    'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.0,
    'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 5.0,
}))

# Procesos (claves de OpenDrift 1.14.9)
for k, v in dict(evaporation=True, emulsification=True, dispersion=True,
                 photooxidation=False, biodegradation=False).items():
    o.set_config(f'processes:{k}', v)
o.set_config('processes:dispersion_mode', a.dispersion_mode)   # NUEVO en v2: legacy | transport | budget
o.set_config('drift:horizontal_diffusivity', 5.0)              # m2/s; en v2 SI se aplica (antes se ignoraba)
o.set_config('general:coastline_action', 'previous')

# Longitud OESTE negativa (lon=+92 cae en el Indico).  max_water=0.9: sin el, MAYA no emulsiona.
o.seed_elements(lon=-92.1, lat=19.7, number=a.n, time=datetime(2023, 7, 6, 12),
                m3_per_hour=100.0, oil_type='MAYA', max_water=0.9)
o.run(duration=timedelta(hours=24), time_step=timedelta(hours=1), outfile=a.out)

b = o.get_oil_budget()
print('DONE', a.out)
print('budget keys:', sorted(b.keys()) if hasattr(b, 'keys') else type(b))
print('Revise el cierre de masa: evaporado + disperso + restante debe sumar la masa inicial.')
