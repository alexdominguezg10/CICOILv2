# CICOILv2

> Traducción al español de `README.md` (la versión en inglés es la referencia).

**Modelo de intemperización de hidrocarburos por componentes para el Golfo de México**

CICOILv2 extiende [OpenDrift](https://github.com/OpenDrift/opendrift) (v1.14.9) con procesos de intemperización de base física resueltos por componente, desarrollados en CICESE (Centro de Investigación Científica y de Educación Superior de Ensenada). Sigue 17 pseudo-componentes a través de evaporación, emulsificación, biodegradación, fotooxidación y disolución submarina, con acoplamiento opcional al modelo de pluma [TAMOC](https://github.com/socolofs/tamoc) para escenarios de descarga submarina.

## Características

- **Balance de masa de 17 componentes** — 8 pseudo-grupos alifáticos + 8 aromáticos + residuo, resueltos a partir de las curvas de destilación de la base de datos ADIOS
- **Biodegradación** — decaimiento de primer orden escalado con Q10, con dependencia del tamaño de gota y vía de fase disuelta
- **Fotooxidación** — degradación de aromáticos inducida por UV con atenuación de Beer-Lambert; seguimiento de fotoproductos oxigenados (OP) con reservorios disuelto/degradado
- **Emulsificación** — captación de agua de Bullwinkle; las constantes opcionales de escalamiento viscosidad-estabilidad (`EMUL_*`) valen 1.0 por defecto (tasa NOAA original)
- **Modos de dispersión** — `processes:dispersion_mode` (`legacy` / `transport` / `budget`) evita contar dos veces la dispersión por rotura de oleaje; sumidero de Roy en la forma de PyGNOME, solo en elementos superficiales
- **Espesor de mancha en malla** — opción `processes:slick_thickness_model=field` (`slick_field.py`): campo de espesor advectado y en esparcimiento que lee la evaporación
- **Acoplamiento con TAMOC** — la salida del modelo de pluma se siembra como elementos lagrangianos
- **Lectores NEMO/CROCO** — lectores NetCDF optimizados para la salida de modelos oceánicos regionales
- **Lector nativo de CROCO** — `reader_croco_native` lee directamente la salida cruda de CROCO en coordenada sigma (`croco_avg.nc` / `croco_his.nc`), con el desescalonamiento correcto de la malla C (véase [Lectura de salida nativa de CROCO](#lectura-de-salida-nativa-de-croco))
- **97 variables de elemento** por partícula

## Requisitos

- Python 3.8+
- OpenDrift 1.14.9
- NumPy, SciPy, xarray, netCDF4, matplotlib
- Base de datos ADIOS (incluida con OpenDrift)
- Opcional: TAMOC (para el acoplamiento con descargas submarinas)

## Instalación

### Instalación para el grupo (Ixachi o laptop Linux): ambiente propio, TAMOC y archivos del modelo que usted puede modificar

Pasos detallados: [`MODEL_CICOILv2/INSTRUCCIONES_INSTALACION_CICOILv2.txt`](MODEL_CICOILv2/INSTRUCCIONES_INSTALACION_CICOILv2.txt). En resumen:

```bash
git clone https://github.com/alexdominguezg10/CICOILv2.git && cd CICOILv2
conda env create -f MODEL_CICOILv2/entorno_cicoilv2.yml        # Python 3.10, numpy 2.x, OpenDrift 1.14.9 (linux-64)
conda activate cicoilv2
bash MODEL_CICOILv2/instalar_cicoilv2.sh preparar               # clona TAMOC (commit f689dc7) en ../tamoc y aplica nuestro parche (necesita git)
bash MODEL_CICOILv2/instalar_cicoilv2.sh instalar               # TAMOC editable + enlaza los archivos de CICOIL de este repositorio en OpenDrift
```

Los archivos de CICOIL se *enlazan* dentro del OpenDrift instalado: al editar `ciceseoil.py` en su clon el cambio se usa de inmediato y aparece en `git diff`.
Se necesitan `gcc` y `gfortran` (TAMOC tiene partes en Fortran); en Ixachi use un nodo de cómputo para el último paso. Verificado desde cero en Ixachi y en una máquina Linux.

### Instalación manual (copiar los archivos en OpenDrift)

CICOILv2 opera como una extensión en sitio de OpenDrift. Instale primero OpenDrift y luego coloque los archivos de CICOILv2 en el directorio del paquete de OpenDrift:

```bash
# 1. Crear el entorno conda
conda create -n cicoil python=3.10
conda activate cicoil

# 2. Instalar OpenDrift 1.14.9
pip install opendrift==1.14.9 future   # 'future' aporta el módulo `past` que usa fluid_properties.py

# 3. Copiar los archivos de CICOILv2 al módulo openoil de OpenDrift
SITE=$(python -c "import opendrift; print(opendrift.__path__[0])")
cp ciceseoil.py fluid_properties.py cicoil_estimations.py \
   tamoc_plume.py run_tamoc.py tamoc_chemical_properties.py \
   cicoil_common.py "$SITE/models/openoil/"
cp -r readers/* "$SITE/readers/"
mkdir -p "$SITE/models/openoil/data"
cp -r data/* "$SITE/models/openoil/data/"
cp export/io_stat_nc.py "$SITE/export/"
```

Para el despliegue en HPC (CHAMAN2), véase `deploy_on_hpc_conda.sh`. En Ixachi existe un entorno compartido; véase §10B del manual. Los scripts listos para correr (Abkatun 2026 en 2D y 3D con TAMOC, Nohoch 2023, verificación de lectores) están en [`examples/ixachi/`](examples/ixachi/).

## Inicio rápido

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

## Lectura de salida nativa de CROCO

`readers/reader_croco_native.py` permite que OpenDrift/CICOILv2 lea **la salida cruda de CROCO** (coordenadas sigma, malla C de Arakawa escalonada) sin preprocesarla a niveles z ni a una malla regular. Es una subclase ligera de `reader_ROMS_native` de OpenDrift 1.14.9 que corrige dos cosas para archivos de CROCO:

1. Los archivos de CROCO traen `mask_rho` pero no `mask_u`/`mask_v` (el lector original falla con `KeyError: 'mask_u'`); ambas se derivan de `mask_rho`.
2. El lector original no desescalona `u` y `v` y los usa como si estuvieran en los nodos rho, un error de posición de media celda (unos 2.3 km en una malla de 1/24 de grado). El nuevo lector promedia `u` y `v` a los puntos rho de forma perezosa (dask) justo después de abrir. `destagger=False` reproduce el comportamiento original.

```python
from opendrift.readers.reader_croco_native import Reader as CrocoNative      # se instala con readers/*
from opendrift.readers.reader_netCDF_CF_generic import Reader as CF

ocean = CrocoNative('croco_avg.nc', name='croco')    # o croco_his.nc; sigma nativo, malla C escalonada
wind  = CF('wind_cf.nc', name='wind')                # los archivos de CROCO no traen atmósfera: dar viento/oleaje aparte
o.add_reader([ocean, wind])
```

Notas:
- No confundir con `reader_NEMO_native_v3.py`, que espera archivos de niveles z estilo NEMO con nombres de variable de NEMO (`vozocrtx`, ...) y no maneja niveles sigma de CROCO.
- Solo variables oceánicas; dé el viento y el oleaje con otro lector.
- El primer registro de `croco_avg.nc` está a la mitad de la primera ventana de promedio (por ejemplo 00:32:30 con promedios horarios), por lo que una corrida no puede empezar antes.
- Validado para **transporte superficial**: frente a una versión convertida de forma independiente a niveles z/malla regular de los mismos archivos, las mismas 5000 partículas se separan 0.2 m en promedio (máx. 2 m) tras 24 h hacia atrás (hindcast horario, trayectoria de ~16 km) y 11 m en promedio (máx. 72 m) tras 48 h hacia adelante (trayectoria de ~87 km); el lector original difiere 2.1 km en promedio (máx. 9.4 km y 3.9 km). Los campos subsuperficiales **aún no están validados**: a 20 m de profundidad las velocidades difieren de esa referencia en una mediana de 0.01-0.04 m/s (causa no establecida).
- Unas cinco veces más lento que leer archivos preprocesados (lee la salida 3-D cruda).
- Escrito y probado con OpenDrift 1.14.9 (se registra una advertencia con otras versiones). Prueba: `python test_croco_native_reader.py` (archivo sintético, 7 verificaciones).

## Estructura de directorios

```
CICOILv2/
├── ciceseoil.py              # Clase principal del modelo (OpenCiceseOil + elemento CiceseOil)
├── fluid_properties.py       # Marco de destilación de 17 componentes
├── cicoil_estimations.py     # Estimación de fracciones SARA
├── tamoc_plume.py            # Pluma TAMOC → siembra en OpenDrift
├── run_tamoc.py              # Envoltura de la simulación TAMOC
├── tamoc_chemical_properties.py
├── cicoil_common.py          # Auxiliares compartidos (load_ciceseoil, default_cicese_config)
├── plotting.py               # Rutinas de graficación para publicación
├── download_data.py          # Descarga y preprocesamiento de GLORYS12/ERA5
├── preprocess_era5_uv.py     # Preprocesamiento de radiación UV de ERA5
├── preprocess_real_forcing.py # Ensamblado de forzamiento de varias fuentes
├── readers/
│   ├── reader_nemo_optimized.py
│   ├── reader_NEMO_native_v3.py
│   ├── reader_nemo_combined.py
│   ├── reader_nemo_modified.py
│   └── reader_croco_native.py      # salida cruda de CROCO en coordenada sigma (desescalonada)
├── export/
│   └── io_stat_nc.py         # Exportación NetCDF estadística en malla
├── data/                     # CSV de propiedades químicas
├── test_cicoil_e2e.py        # Batería de extremo a extremo de 30 pruebas (T1–T8)
├── test_croco_native_reader.py  # reader_croco_native en un archivo CROCO sintético (7 verificaciones)
├── test_horizontal_diffusivity.py  # paseo aleatorio de drift:horizontal_diffusivity (4 verificaciones)
├── test_biodegradation.py
├── test_biodegradation_water.py
├── test_droplet_size_biodeg.py
├── test_emulsification_viscosity.py
├── test_op_tracking.py
├── test_photoox_era5uv.py
├── test_photoox_vertmix.py
├── deploy_on_hpc_conda.sh    # Script de despliegue en HPC
├── CICOIL_v2_MANUAL.md       # Manual completo (física + despliegue)
├── CICOIL_v2_MANUAL_ES.md    # Manual completo en español
├── CICOIL_v2_REFERENCE.md    # Referencia de la API
└── BIODEGRADATION.md         # Documentación del módulo de biodegradación
```

## Pruebas

```bash
conda run -n cicoil python test_cicoil_e2e.py
# Esperado: PASSED: 30/30
```

La batería cubre:
- T1–T3: modos de intemperización NOAA y CICESE, conservación de masa
- T4: siembra de pluma TAMOC
- T5: importaciones de lectores NEMO
- T6: graficación del balance de masa
- T7: Biodegradación (balance de masa por componente)
- T8: Fotooxidación (seguimiento de OP + cotas de regresión)

## Documentación

- **[CICOIL_v2_MANUAL_ES.md](CICOIL_v2_MANUAL_ES.md)** — Manual completo: instalación, despliegue, física de la intemperización, calibración de parámetros
- **[CICOIL_v2_REFERENCE_ES.md](CICOIL_v2_REFERENCE_ES.md)** — Referencia de la API: variables de elemento, claves de configuración, firmas de métodos
- **¿Viene de CICOIL v1?** Lea la §10B del manual (qué cambia, por qué difieren los resultados, instalación compartida en Ixachi). La §10C lista los cambios posteriores a `b032140` (máscara de burbujas de gas, corrección de ID en `oil_density`, `vertical_profile.py`).

## Autores

- **Konstantinos Kotzakoulakis** (CICESE / SINTEF Ocean) — desarrollo original de CICOILv2
- **Knut-Frode Dagestad** (MET Norway) — marco de OpenDrift
- **Alex Dominguez** (CICESE) — biodegradación, fotooxidación, seguimiento de OP, correcciones de código

## Licencia

GNU General Public License v2.0 — véase [LICENSE](LICENSE).

CICOILv2 es una obra derivada de [OpenDrift](https://github.com/OpenDrift/opendrift), distribuida bajo la misma licencia.
