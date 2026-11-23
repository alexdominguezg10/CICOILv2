# Manual de CICOILv2

> Traducción al español de `CICOIL_v2_MANUAL.md` (versión en inglés = referencia; si difieren, manda el inglés).

**Modelo de intemperización de hidrocarburos basado en OpenDrift (fork CICOIL, OpenDrift 1.14.9 + extensiones de CICESE)**

Este manual reúne y sustituye a `cicoilv2/hpc_migration.md` (flujo de
despliegue) y `DEPLOYMENT_MANIFEST.txt` (inventario de archivos / referencia de
la API), y agrega los módulos de fotooxidación, biodegradación y emulsificación
desarrollados hasta el 2026-06-15 y los cambios de intemperización de octubre de
2026 (§8A, §8B, §10A).

**Copia de desarrollo activa:** `OILSPILL/CICOIL_dev/cicoilv2_deployment/`
**Contacto:** adomingu@cicese.edu.mx — CICESE, Departamento de Oceanografía Física

> Para un catálogo completo de la API y la funcionalidad (cada clase, método,
> función, variable de elemento y clave de configuración, incluidos FluidProps,
> el acoplamiento con TAMOC, los lectores NEMO/CROCO, y las utilidades de
> graficación y descarga), véase el documento complementario
> [`CICOIL_v2_REFERENCE_ES.md`](CICOIL_v2_REFERENCE_ES.md).

---

## 1. Resumen

CICOILv2 = `OpenOil` de OpenDrift 1.14.9 + las extensiones de intemperización de
CICOIL (Kotzakoulakis, CICESE/SINTEF Ocean), parchadas para ser compatibles con
OpenDrift 1.14.9, más tres mejoras desarrolladas en CICESE:

- **Fotooxidación** (`photooxidation_cicese()`) — pérdida, por radiación UV
  solar, de los pseudo-componentes aromáticos en la superficie o cerca de ella.
- **Biodegradación** (`biodegradation_cicese()`) — resuelta por componente,
  escalada con la temperatura, con supresión en la mancha superficial y una vía
  de fase disuelta.
- **Emulsificación** (`emulsification_cicese()`) — captación de agua de
  Bullwinkle con escalamiento opcional por viscosidad-estabilidad (constantes
  `EMUL_*`). La calibración de bolas de alquitrán de junio de 2026 se retiró el
  2026-10-29: las constantes valen 1.0 (tasa NOAA original); véase §8.

La fotooxidación y la biodegradación son aditivas (apagadas por defecto) y
conservan la masa exactamente (error relativo ~1e-15–1e-16 en todas las
pruebas). La emulsificación está activa con `processes:emulsification = True`.
Los cambios desde v2.0.0 están en §10A; quienes vienen de CICOIL v1 deben leer
§10B.

### Qué agrega CICOIL a `OpenOil` de OpenDrift 1.14.9
- Intemperización de 17 componentes (evaporación por presión de vapor de Cox,
  esparcimiento de Nihoul, densidad por partícula)
- Acoplamiento con el modelo de pluma TAMOC (`seed_plume_elements()`, `Plume`,
  `run_tamoc()`)
- Modelo de dispersantes (superficial + SSDI, 4 tipos de dispersante, IFT y
  max_water por partícula)
- 3 tipos de partícula: gotas superficiales, burbujas de gas (`is_bubble`) y
  disuelto submarino (`is_dissolved`)
- 93 variables de elemento (frente a 22 en `OpenOil` original de OD 1.14.9)
- 4 lectores NEMO/CROCO (native_v3, combined, modified, optimized)
- Salida NetCDF estadística (`export/io_stat_nc.py`)
- **NUEVO (CICESE):** fotooxidación, biodegradación y emulsificación,
  descritas en §6–8

---

## 2. Instalación y despliegue

### 2.1 Entorno de desarrollo local
- Árbol de fuentes: `OILSPILL/CICOIL_dev/cicoilv2_deployment/`
- `OILSPILL/CICOIL` es un enlace simbólico a la ubicación de la instalación de
  OpenDrift
  (`~/miniconda3/envs/numerics/lib/python3.11/site-packages/opendrift/`)
- Entorno conda local: `numerics` (Python 3.11, OpenDrift 1.14.9, TAMOC 4.x)

Para correr las pruebas localmente:
```bash
conda run -n numerics python test_cicoil_e2e.py
conda run -n numerics python test_biodegradation.py
conda run -n numerics python test_biodegradation_water.py
```

### 2.2 Despliegue en HPC (CHAMAN2)

**Entorno destino:** `cicoil4` (Python 3.10, clonado de `cicoil3`).
**Gestor de paquetes:** conda (mamba no funciona en CHAMAN2: no use
`deploy_on_hpc_mamba.sh`, use **`deploy_on_hpc_conda.sh`**).
**Seguridad:** `cicoil3` nunca se modifica; `cicoil4` es un entorno nuevo y
dedicado que se puede borrar y recrear.

#### Paso 1 — Preparación local
```bash
cd OILSPILL/CICOIL_dev/cicoilv2_deployment
bash stage_for_hpc.sh   # genera cicoilv2_hpc_deployment.tar.gz + DEPLOYMENT_MANIFEST.txt
```

#### Paso 2 — Transferencia
```bash
scp cicoilv2_hpc_deployment.tar.gz YOUR_USERNAME@chaman2.cicese.mx:~/cicoilv2_migration/
ssh YOUR_USERNAME@chaman2.cicese.mx
cd ~/cicoilv2_migration
tar -xzf cicoilv2_hpc_deployment.tar.gz
cd cicoilv2_deployment
```

#### Paso 3 — Despliegue
```bash
bash deploy_on_hpc_conda.sh
```
Este script:
1. Crea `cicoil4` (clonado de `cicoil3`, conda)
2. Actualiza OpenDrift 1.9.0 → 1.14.9 solo en `cicoil4`
3. Instala todos los archivos parchados/nuevos en el site-packages de `cicoil4`
4. Corre la batería completa de pruebas (`test_cicoil_e2e.py`,
   `test_biodegradation.py`, `test_biodegradation_water.py`)
5. Escribe `HPC_DEPLOYMENT_REPORT.txt`

Resumen esperado al final del script:
```
=== HPC DEPLOYMENT REPORT ===
✓ OpenDrift version: 1.14.9
✓ Patched files installed
✓ test_cicoil_e2e.py: 30/30 PASSED
✓ test_biodegradation.py: PASSED (Q10 sanity check)
✓ test_biodegradation_water.py: PASSED
✓ Mass conservation: 100.000% (rel. err ~1e-15)
✓ Ready for production — activate with: conda activate cicoil4
```

#### Paso 4 — Corrida de producción
```bash
cp slurm_template_cicoilv2.sh ~/my_simulation.slurm
nano ~/my_simulation.slurm   # poner conda activate cicoil4, editar los datos del experimento
sbatch ~/my_simulation.slurm
squeue -u YOUR_USERNAME
```

#### Paso 5 — Verificar la salida
```python
import xarray as xr
ds = xr.open_dataset('output.nc')
print(len(ds.time), ds.sizes['trajectory'])
print(set(ds.status.values.ravel()))
```

---

## 3. Inicio rápido

```python
import logging; logging.disable(logging.WARNING)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_netCDF_CF_generic import Reader as CFReader
from datetime import timedelta

o = OpenCiceseOil(weathering_model='cicese')
o.set_config('processes:evaporation', True)
o.set_config('processes:spreading', True)
o.set_config('processes:emulsification', True)
o.set_config('processes:photooxidation', True)
o.set_config('processes:biodegradation', True)

o.add_reader([CFReader('glorys12_ocean.nc'), CFReader('era5_wind_uv_cf.nc')])

o.seed_elements(lon=-92.1, lat=19.7, time='2023-07-06 12:00',
                number=3000, m3_per_hour=100.0,
                oil_type='MAYA', max_water=0.9)

o.run(duration=timedelta(days=30), time_step=timedelta(hours=1),
      time_step_output=timedelta(hours=6), outfile='output.nc')

o.plot_oil_budget(filename='oil_budget.png')
```

**Tamaño de la salida en ensambles grandes:** por defecto `o.run()` omite de
`output.nc` los 51 arreglos por pseudo-componente (`comp_*`,
`comp_*_evaporated`, `comp_*_dissolved`), que dominan el tamaño del archivo
(17 componentes x 3 variantes x N partículas x N pasos de salida). Con
`record_components=True` se conservan (p. ej. para depurar el balance de masa).
En ensambles muy grandes, pase `record_environment=False` para omitir también
los 20 campos de ambiente por elemento (corriente, viento, oleaje, SST,
salinidad, etc.); se pueden rederivar de los datos de forzamiento y suelen ser
la mayor contribución restante al tamaño. OpenDrift comprime la salida NetCDF
con zlib (complevel 6) al cerrar el archivo, con o sin estas opciones.

**Descarga submarina TAMOC (siembra de pluma):**
```python
from opendrift.models.openoil.tamoc_plume import Plume
from opendrift.models.openoil.run_tamoc import run_tamoc
# construir Plume a partir de la salida de run_tamoc(...), y luego:
o.seed_plume_elements(plume, ...)
```

**Forzamiento CROCO/NEMO:**
```python
from opendrift.readers.reader_nemo_optimized import Reader as CrocoReader
o.add_reader(CrocoReader('croco_avg.nc'))
```

**Salida cruda de CROCO en coordenada sigma (sin preprocesamiento):**
```python
from opendrift.readers.reader_croco_native import Reader as CrocoNative      # se instala con readers/*
from opendrift.readers.reader_netCDF_CF_generic import Reader as CF

ocean = CrocoNative('croco_avg.nc', name='croco')    # o croco_his.nc; sigma nativo, malla C escalonada
wind  = CF('wind_cf.nc', name='wind')                # los archivos de CROCO no traen atmósfera: dar viento/oleaje aparte
o.add_reader([ocean, wind])
```
`reader_croco_native` hereda de `reader_ROMS_native` de OpenDrift, deriva
`mask_u`/`mask_v` de `mask_rho` y desescalona `u`/`v` a los puntos rho
(`destagger=True`, por defecto; `False` = comportamiento original con un
desplazamiento de media celda). Lee `croco_avg.nc` / `croco_his.nc` tal como
salen de CROCO; a diferencia de `reader_NEMO_native_v3`, maneja coordenadas
sigma y la malla C. Solo variables oceánicas; el primer registro de
`croco_avg.nc` está a la mitad de la primera ventana de promedio, por lo que una
corrida no puede empezar antes. Validado para transporte superficial (24 h hacia
atrás / 48 h hacia adelante, 5000 partículas: separación media de 0.2 m / 11 m
frente a archivos convertidos, contra 2.1 km del lector original); los campos
subsuperficiales aún no están validados (diferencias de velocidad de 0.01-0.04
m/s a 20 m de profundidad contra la referencia). Probado con OpenDrift 1.14.9.

### Verificación rápida
```python
from opendrift.models.openoil.ciceseoil import OpenCiceseOil, CiceseOil
o = OpenCiceseOil(weathering_model='cicese')
o.set_oiltype('MAYA')
assert len(CiceseOil.variables) == 93
assert o.get_config('biodegradation:method') == 'Adcroft'
print("CICOIL OK:", len(o.oiltypes), "oils,", len(o.fluid_properties.composition), "components")
```

---

## 4. Referencia de configuración (claves nuevas o modificadas)

| Clave | Valor por defecto | Notas |
|---|---|---|
| `processes:evaporation` | True | Evaporación por Pv de Cox |
| `processes:spreading` | True | Esparcimiento de Nihoul |
| `processes:emulsification` | True | Modelo de Bullwinkle (escalamiento viscosidad-estabilidad disponible con `EMUL_*`, neutro con el valor 1.0; véase §8); requiere `max_water>0` al sembrar (véase la advertencia de §9) |
| `processes:dispersion` | False | Sumidero de masa por dispersión natural de Roy (solo elementos superficiales; véase `processes:dispersion_mode`, §8A) |
| `processes:dispersion_mode` | `legacy` | NUEVO — `legacy` / `transport` / `budget`: cuál de los dos, el sumidero de Roy o el arrastre de Li, actúa (§8A) |
| `processes:slick_thickness_model` | `legacy` | NUEVO — `legacy` / `field`: espesor de mancha en malla que lee la evaporación (§8B); `field` requiere `processes:spreading` False |
| `slick_field:*` | véase Reference §3.3 | NUEVO — caja de la malla (por defecto Golfo de México), `dx` 500 m, `lens_age`, `evaporation_thickness`, `drift_update_interval` |
| `processes:chemical_dispersion` | False | Modelo de dispersantes superficial + SSDI |
| `processes:subsea_dissolution` | True | Alimenta el reservorio de fase disuelta (§7.3) |
| `processes:surface_dissolution` | False | |
| `processes:update_oilfilm_thickness` | True | |
| `processes:handle_released_gas` | True | |
| `processes:photooxidation` | **False** | NUEVO — §6 |
| `processes:biodegradation` | **False** | NUEVO — §7 |
| `biodegradation:method` | `'Adcroft'` | Ahora describe la ley Q10 por componente de `biodegradation_cicese()`, no la fórmula global original |
| `seed:oil_type` | — | Nombre del crudo en ADIOS, p. ej. `'MAYA'` |
| `seed:m3_per_hour` | — | Requerido para siembra continua |
| `seed:droplet_size_distribution` | — | `'normal'`/`'lognormal'`/`'uniform'` para liberación submarina |
| `seed:droplet_diameter_mu` / `_sigma` | — | Parámetros de la distribución de tamaños de gota submarina |
| `seed:droplet_diameter_min/max_subsea` | — | |
| `drift:profiles_depth` | 20 | Nombre de clave en plural (OD 1.14.9) |
| `drift:vertical_advection` | False | |
| `drift:vertical_mixing` | True | Necesario para el comportamiento correcto de la fracción superficial en fotooxidación (§6.2) |
| `drift:horizontal_diffusivity` | 0 | Difusividad del paseo aleatorio [m² s⁻¹]. En OpenDrift 1.14.9 esta clave nunca se lee y `horizontal_diffusion()` se omite a menos que `horizontal_diffusivity` esté en `required_variables`; CICOIL agrega esa variable (valor de respaldo 0) y refleja esta clave en `environment:fallback:horizontal_diffusivity`. Dispersión ≈ √(2 D t). Las corridas hechas antes de esta corrección no tenían difusión, fuera cual fuese el valor. |

---

## 5. Correcciones de compatibilidad de la API aplicadas (14, frente a OD 1.14.9)

1. `adios.computation` → `adios_db.computation` (importaciones)
2. `adios.util.estimations` → `cicoil_estimations.py` independiente
3. `adios.oil_name_alias` → protegido con `getattr()`
4. `CONFIG_LEVEL_*` ahora se importa de `opendrift.config` (a nivel de módulo, no atributos de clase)
5. atributo de clase `max_water_fraction` → `_cicoil_max_water_fraction_override`; `self.max_water_fraction = None` en `__init__`
6. `self.get_environment()` → `self.env.get_environment()` (devuelve `(recarray, profiles, missing)`)
7. `self.Density` / `self.KinematicViscosity` se inicializan en `seed_elements()` en modo NOAA
8. Corregida la indexación del ID de elemento: `elements.ID - 1` → `elements.ID` (4 ocurrencias)
9. `__set_seed_config__()` para `seed:oil_type` (funciona en los modos Config y Ready; corrige la guarda de modo de `seed_plume_elements()`)
10. Corregida la transición Mode.Config → Mode.Ready en `seed_plume_elements()`
11. `pyproj` se importa directamente (no vía `basereader.pyproj`)
12. `get_oil_budget()` corregido para xarray/arreglos enmascarados mediante la función auxiliar `_prop()` (`.values` + `np.nansum()`)
13. se quitó `required_profiles_z_range`; se reemplazó por `_set_config_default('drift:profiles_depth', 20)`
14. respaldo de finalización del ambiente en `seed_plume_elements()` (temperatura por defecto de 28 °C antes de finalizar el ambiente)

---

## 6. Módulo de fotooxidación (`photooxidation_cicese()`)

Los grupos aromáticos de pseudo-componentes G1–G8 se fotooxidan por la
radiación UV solar cerca de la superficie. **Estado: CERRADO el 2026-06-13**:
se considera suficiente; no se planea más desarrollo salvo que haya nuevos
datos de validación.

### 6.1 Modelo de tasa
```
dM_i = K_PHOTO_AROM[i] * I_UV(z) * dt * M_i      (i = G1..G8, solo aromáticos)
I_UV(z) = I_UV(0) * exp(-K_D_UV * |z|)            (Beer-Lambert, resuelto en profundidad)
```
- Los elementos más profundos que `UV_DEPTH_CUTOFF = 50 m` se omiten.
- `K_D_UV = 0.10 m⁻¹` (condición de "Golfo de México claro", de `surface_fate_coupled.py`).
- Componentes alifáticos y residuo: no se afectan (0.000 kg).
- Solo actúa de día (`I_UV(0) = 0` de noche).

### 6.2 Forzamiento UV
- **Preferido**: `surface_downward_uv_radiation` real de ERA5 (preprocesado a CF
  con `preprocess_era5_uv.py`, datos en `data_local/era5_uv_cf/`, 1972–2026).
  El UV real de ERA5 llega a ~120 W/m² (frente a ~50 W/m² del sintético anterior).
- **Respaldo**: modelo sintético de coseno diurno del ángulo cenital solar (pico
  de 50 W/m², Bahía de Campeche; solo si no se da un lector de UV).
- Requiere `drift:vertical_mixing = True`: la mezcla vertical 3D de OpenDrift ya
  representa la fracción superficial; **no** aplique además una corrección
  `f_surface` aparte (la antigua `_compute_f_surface_cicese()` se eliminó por
  contar doble). Con `vertical_mixing=False` la fotooxidación se sobrestima ~11 %.

### 6.3 Calibración actual (`K_PHOTO_AROM`, m² J⁻¹)
| Grupo | Compuesto | k (m² J⁻¹) | t½ |
|---|---|---|---|
| G1 | benceno | 2.9e-8 | ~14 d |
| G2 | tolueno | 2.0e-8 | ~20 d |
| G3 | xilenos | 5.0e-8 | ~8 d |
| G4 | naftaleno | 4.260e-8 | ~9.4 d |
| G5 | acenafteno | 5.145e-8 | ~7.8 d (interpolado G4–G6) |
| G6 | fenantreno | 6.176e-8 | ~6.5 d |
| G7 | criseno | 1.581e-8 | ~25.3 d (interpolado G6–G8) |
| G8 | pireno+ | 5.110e-9 | ~79.8 d |

G4/G6/G8 se recalibraron por mínimos cuadrados de la serie completa (f\*=0.3834)
contra la referencia a 30 días de `surface_fate_coupled.py`
(`data_local/surface_fate/surface_fate_base.csv`). G5/G7 son interpolaciones
log-lineales contra el punto de ebullición (no se encontró literatura específica
por compuesto; los valores caen dentro de los rangos por clase de anillo de
Bacosa et al. 2015). Persiste una discrepancia de forma: el esquema resuelto en
profundidad se satura hacia los días 10–15, mientras que el objetivo sube
gradualmente en 30 días (limitación documentada, no corregida).

### 6.4 Resultado de referencia e incertidumbre
Nohoch Alfa 2023 (N=3000, 30 d, forzamiento real, MAYA): **total_oxy_frac ≈
5.73 %** (~31.9 t de un derrame de 557 t). Envolvente de incertidumbre
(un parámetro a la vez: K_PHOTO_AROM ±50 %, K_D_UV 0.05–0.20, DELTA_M_OX 16–48
g/mol): **2.85–9.60 %** (~15.9–53.5 t, 0.5×–1.7× el valor base): cota
conservadora, no estadística. K_PHOTO_AROM domina la incertidumbre (factor ~2.4).

**Interpretación**: 5.73 % es una estimación de **cota inferior**: solo
aromáticos (G1–G8), sin vía de fotooxidación alifática (cf. Ward et al. 2018
sobre fotooxidación alifática mediada por radicales).

### 6.5 Variables de elemento nuevas
`mass_photooxidized` (kg), `fraction_photooxidized` (%). `mass_photooxidized` es
un **sumidero terminal**: no se biodegrada más.

### 6.6 Seguimiento de fotoproductos oxigenados (OP) (auditoría #4, 2026-06-14)
De la masa fotooxidada en cada paso `dM_phox` (grupos aromáticos G1–G8), una
fracción se dirige a un reservorio disuelto de fotoproductos oxigenados en vez
del sumidero de fase oleosa `mass_photooxidized`:
```
dM_op      = OP_DISSOLVED_FRACTION * dM_phox        # -> mass_op_dissolved
dM_residue = (1 - OP_DISSOLVED_FRACTION) * dM_phox  # -> mass_photooxidized (como antes)
```
- `OP_DISSOLVED_FRACTION = 0.5`: reparto parejo 50/50, estimación de primera
  aproximación que refleja que los OP (oxi-HAP, quinonas, ácidos carboxílicos)
  son notablemente más solubles en agua que sus HAP padres (Aeppli et al. 2012);
  pendiente de datos de partición por clase de compuesto.
- `mass_op_dissolved` decae con su propia tasa de primer orden
  `K_OP_REMOVAL = 0.116 d⁻¹` (= 2 × `K_BIODEG_WATER_ANCHOR`, t½≈6 d a 20 °C,
  escalada con Q10 como `K_BIODEG_*`), hacia `mass_op_degraded`. Este paso de
  remoción es **incondicional** (corre incluso de noche / `I_UV(0)=0`), porque
  representa un proceso de la columna de agua (biodegradación, fotodegradación
  o dilución posteriores de los OP) independiente del UV solar.
- `fraction_photooxidized`, `fraction_op_dissolved`, `fraction_op_degraded`
  comparten el mismo denominador conservado `initial_mass` (que ahora incluye
  también `mass_op_dissolved` + `mass_op_degraded`); el denominador análogo en
  `biodegradation_cicese()` (`fraction_biodegraded_from_oil`) se actualizó igual
  por consistencia.
- Variables de elemento nuevas: `mass_op_dissolved`, `fraction_op_dissolved`,
  `mass_op_degraded`, `fraction_op_degraded` (todas kg/%, por defecto 0).
- Verificado por `test_op_tracking.py`: corrida superficial de MAYA de 8 días,
  conservación de masa incluidos los nuevos reservorios (error relativo 0),
  `mass_op_dissolved + mass_op_degraded == mass_photooxidized` (exacto, por la
  construcción del reparto 50/50), `mass_op_degraded > 0` (el decaimiento
  K_OP_REMOVAL está activo).

---

## 7. Módulo de biodegradación (`biodegradation_cicese()`)

Agregado el 2026-06-14, en lugar de la vía anterior, inactiva e inconsistente
(el código previo llamaba a `biodegradation()` heredado →
`biodegradation_adcroft()`, que hacía decaer `mass_oil` global sin actualizar
`cicese_mass_balance['mass_components']`: una inconsistencia del balance de masa;
las variables de seguimiento `mass_biodegraded_from_oil/_water` de CICOIL
existían pero siempre valían cero).

### 7.1 Modelo de tasa
Decaimiento exponencial de primer orden resuelto por componente, aplicado a cada
pseudo-componente:
```
dM = M * (1 - exp(-K(T) * dt))
K(T) = K(T_REF) * BIODEG_Q10 ** ((T - T_REF) / 10)     # T_REF=20°C, BIODEG_Q10=3
```
- Se usa `sea_water_temperature` por elemento (con soporte de perfiles:
  funciona para gotas submarinas).
- Actúa sobre **todas las partículas activas** (sin límite de profundidad ni de
  luz, a diferencia de la fotooxidación).
- Las partículas superficiales (`z==0`) usan `K * BIODEG_SURFACE_FACTOR` con
  **`BIODEG_SURFACE_FACTOR = 1/30`**: Prince et al. (2003, 2017) encontraron
  que las manchas flotantes están "casi inmunes a la biodegradación detectable"
  frente al crudo disperso (vida media de 1–3 semanas), con una vida media
  superficial aparente de "muchos meses a muchos años". `1/30` lleva, p. ej., una
  vida media de 5 d en crudo disperso a ~5 meses en superficie.

### 7.2 Calibración actual (`K_BIODEG_*`, d⁻¹ a T_REF=20 °C, crudo disperso/arrastrado z<0)
| Grupo | k alifático | t½ | k aromático | t½ |
|---|---|---|---|---|
| G1 | 0.1386 | 5 d | 0.0866 | 8 d |
| G2 | 0.1155 | 6 d | 0.0693 | 10 d |
| G3 | 0.0990 | 7 d | 0.0578 | 12 d |
| G4 | 0.0866 | 8 d | 0.0462 | 15 d |
| G5 | 0.0693 | 10 d | 0.0347 | 20 d |
| G6 | 0.0462 | 15 d | 0.0277 | 25 d |
| G7 | 0.0277 | 25 d | 0.0173 | 40 d |
| G8 | 0.0173 | 40 d | 0.0116 | 60 d |

Residuo: `K_BIODEG_RESIDUE = 0.00139` (t½ ~500 d, prácticamente no biodegradable
en las escalas de tiempo de simulación). Partículas superficiales: multiplicar
todo lo anterior por `BIODEG_SURFACE_FACTOR = 1/30`.

Base: ley Q10 de Adcroft et al. (2010); magnitudes relativas de Bacosa et al.
(2015) (los n-alcanos ~10× más rápidos que la fotooxidación de HAP; biodegradación
comparable o más rápida que la fotooxidación para HAP pesados de 4–5 anillos);
vidas medias alifáticas de mesocosmos de 8–19 d. **De primera aproximación,
informada por la literatura; pendiente de recalibración con datos** (mismo
estado que K_PHOTO_AROM antes de recalibrar).

### 7.3 Vía de fase disuelta (`mass_biodegraded_from_water`)
`handle_subsea_dissolution()` ahora mantiene un reservorio persistente de masa
disuelta, en lugar de desactivar de inmediato:
- En el primer encuentro (`is_dissolved==1`): transferencia única de
  `cicese_mass_balance['mass_components']` → `mass_dissolved` (componentes en
  cero); `elements.mass_dissolved_subsea` se sobrescribe con este total de
  referencia.
- En cada paso: `mass_dissolved` decae a **`K_BIODEG_WATER = 0.0578 d⁻¹`**
  (t½~12 d a 20 °C, escalada con Q10, **sin** supresión superficial: tasa de
  referencia de pluma disuelta del Golfo de México de Adcroft et al. 2010,
  también consistente con el extremo rápido del rango de 1–3 semanas de Prince
  et al. para crudo disperso).
- El elemento se desactiva (`reason='dissolved_biodegraded'`) cuando
  `mass_dissolved < BIODEG_WATER_DEPLETION_FRAC (=0.01, es decir 1 %)` de la
  masa de referencia inicial.
- Tasa única de primera aproximación, no resuelta por componente.
- Solo relevante para escenarios de descarga submarina con TAMOC
  (`seed_plume_elements()`); no se ejercita en liberaciones superficiales
  (`seed_elements()`).

### 7.4 Variables de elemento nuevas/actualizadas
`mass_biodegraded_from_oil`, `mass_biodegraded_from_water`, el genérico
`mass_biodegraded` (para reportes/gráficas existentes),
`fraction_biodegraded_from_oil`, `fraction_biodegraded_from_water`. Actualiza
el `cicese_mass_balance['mass_components']` por componente (mantiene el balance
por componente consistente con `mass_oil`).

### 7.5 Resultados de referencia — Nohoch Alfa 2023 (N=3000, 30 d, forzamiento real, MAYA, liberación superficial)
% de la masa inicial total de crudo al día 30, con `processes:biodegradation` y
`processes:photooxidation` activados:

| Vía | Antes de recalibrar (sin supresión superficial) | Recalibrado (BIODEG_SURFACE_FACTOR=1/30) |
|---|---|---|
| Remanente | 3.0329 % | 4.4324 % |
| Evaporado | 73.8170 % | 81.5155 % |
| Fotooxidado | 3.2750 % | 4.2874 % |
| **Biodegradado (del crudo)** | **19.8751 %** | **9.7648 %** |
| Error relativo de conservación de masa | 0.00e+00 | 1.36e-16 |

Escalado al tamaño del derrame: 557 t → 54.4 t biodegradadas / 23.9 t
fotooxidadas; 1000 t → 97.6 t / 42.9 t.

### 7.6 Biodegradación resuelta por tamaño de gota (auditoría #2, 2026-06-14)
La matriz de tasas por grupo `K` (después del escalamiento por Q10 y supresión
superficial) se multiplica además por un factor de tamaño de razón
superficie/volumen:
```
K_i_eff = K_i * (BIODEG_DROPLET_D_REF / d) ** BIODEG_DROPLET_EXPONENT
```
- `BIODEG_DROPLET_D_REF = 0.001 m` (1 mm): coincide con
  `seed:droplet_diameter_mu`, tomado como diámetro de referencia de calibración
  de `K_BIODEG_ALI/AROM/RESIDUE` (§7.2).
- `BIODEG_DROPLET_EXPONENT = 1.0`: SA/V ~ 1/d para una esfera (King 1992;
  Brakstad et al. 2015): las gotas más pequeñas presentan proporcionalmente más
  interfaz aceite-agua por unidad de volumen de crudo y se degradan más rápido.
- `d = elements.diameter`, acotado por `BIODEG_DROPLET_DIAMETER_MIN/MAX =
  [1e-5, 1e-2] m` (10 µm – 1 cm), lo que acota `K_i_eff/K_i` a `[0.1, 100]`.
- `diameter == 0` (el valor por defecto de la mayoría de las partículas
  superficiales "spillet", para las que `diameter` no es un tamaño de gota con
  significado) se trata como `d_ref`, es decir factor = 1: sin cambio y sin
  contar doble con `BIODEG_SURFACE_FACTOR` (§7.1), que sigue siendo el
  discriminador principal superficie-disperso.

Efecto neto: las gotas submarinas/dispersas menores de 1 mm (p. ej. arrastre
fino por oleaje o gotas de descarga submarina) se biodegradan más rápido que la
línea base de §7.2; las gotas mayores de 1 mm, más lento. Las manchas
superficiales (`diameter==0`) no se afectan.

Verificado por `test_droplet_size_biodeg.py` (90 partículas, 4 grupos de 30:
`diameter` = 0, 0.1 mm, 1 mm, 1 cm, 10 d a SST=20 °C=`BIODEG_T_REF`): error
relativo de conservación de masa < 1e-6;
`frac_biodeg(0.1mm)=7.110% > frac_biodeg(1mm)=0.810% > frac_biodeg(1cm)=0.082%`;
`frac_biodeg(diameter=0) == frac_biodeg(diameter=1mm)` (dif. relativa 4.7e-8).

---

## 8. Módulo de emulsificación (`emulsification_cicese()`)

### 8.1 Problema con la tasa de Bullwinkle original de NOAA

El modelo original de captación de agua de Bullwinkle de NOAA PyGNOME calcula
una única escala de tiempo característica `t_char = (6/drop_max) / k_emul_eff`
donde, con un viento típico de 10 m/s, el `k_emul ≈ 121 s⁻¹` sin escalar da
`t_char ≈ 1–5 horas`: el crudo forma una mousse estable en pocas horas sin
importar su viscosidad ni su estado de intemperización. Para el crudo MAYA (y la
mayoría de los crudos medios a ligeros), las observaciones de campo y Fingas &
Fieldhouse (2004) indican que la formación de bolas de alquitrán/mousse toma
**días a semanas**, no horas.

### 8.2 Escalamiento viscosidad-estabilidad

La corrección introduce un factor de estabilidad dependiente de la viscosidad
antes de aplicar la tasa de captación de agua de Bullwinkle de NOAA:

```
k_emul_eff = k_emul * clip(nu_dry / EMUL_STABILITY_NU_REF,
                            EMUL_STABILITY_MIN, EMUL_STABILITY_MAX)
           * EMUL_RATE_SCALE
```

donde:
- `nu_dry` = viscosidad cinemática [cSt] del crudo seco (sin emulsión), corregida
  por evaporación (`self._oil_viscosity_dry`, calculada en cada paso en
  `oil_weathering_cicese()`).
- `EMUL_STABILITY_NU_REF = 1000.0 cSt`: viscosidad de referencia en torno a la
  cual se centra el factor de estabilidad (el argumento de clip vale 1.0 con
  `nu_dry = 1000 cSt`).
- `EMUL_STABILITY_MIN = 1.0`: factor de estabilidad mínimo (era 0.3).
- `EMUL_STABILITY_MAX = 1.0`: factor de estabilidad máximo (era 3.0).
- `EMUL_RATE_SCALE = 1.0`: corrección multiplicativa global (era 0.0142).

> **Actualización 2026-10-29.** Las tres constantes vuelven a valer 1.0
> (atributos de clase y `CICESE_DEFAULTS`), de modo que el valor por defecto es
> la tasa de Bullwinkle original de NOAA, independiente de la viscosidad. El
> ajuste descrito en §8.2–§8.4 (0.3 / 3.0 / 0.0142, ajustado a un objetivo de
> bolas de alquitrán a 17.2 días para Abkatun 2026) no tenía respaldo en datos:
> la captación de agua era unas 70 veces más lenta que la tasa original
> (fracción de agua de 0.12 a las 48 h en las corridas S0 de Nohoch 2023, contra
> 0.78 a las 12 h en v1), y la viscosidad se mantenía cerca de 1e-3 m² s⁻¹, de
> modo que el sumidero de dispersión de Roy nunca se apagaba (modo budget: 22 %
> disperso contra 4.4 % en PyGNOME). Los resultados citados en §8.3–§8.4 y en §9
> se obtuvieron con los valores anteriores; las corridas que fijan las
> constantes en tiempo de ejecución no se afectan. Para reproducirlos, fije
> `EMUL_STABILITY_MIN = 0.3`, `EMUL_STABILITY_MAX = 3.0`,
> `EMUL_RATE_SCALE = 0.0142` en la clase antes de sembrar.

### 8.3 Calibración empírica (histórica; constantes revertidas el 2026-10-29, véase la actualización anterior)

Calibrada contra una corrida realista TAMOC+GLORYS12 de Abkatun 2026 con crudo
MAYA (`calibrate_emulsification_viscosity_tamoc.py`), usando un objetivo de
transición a bolas de alquitrán superficiales a 17.2 días de
`TarballFormation_MPB.tex` (consistente con Fingas & Fieldhouse 2004):

| `EMUL_RATE_SCALE` | Día de cruce a bola de alquitrán (ν ≥ 10,000 cSt) |
|---|---|
| 1.0 (sin escalamiento) | 0.76 d |
| 0.024 | 10.21 d |
| **0.0142** | **16.92 d** ← elegido (objetivo 17.2 d, error <1.6 %) |
| 0.0072 | nunca en 25 d (ν = 7,388 cSt al día 25) |

### 8.4 Impacto en los balances de masa de los casos de estudio

La corrección de emulsificación fue el **cambio dominante** en las re-corridas
de Nohoch 2023 y Abkatun 2026 (junio de 2026):

- **Antes de la corrección**: toda partícula alcanzaba `max_water=0.9` (el tope
  de ADIOS) y una viscosidad 18–32× por encima del umbral de bola de alquitrán
  en unas ~5 horas desde la liberación, en todos los casos.
- **Después de la corrección**: el contenido de agua aumenta gradualmente en
  días, llegando a `water_fraction` ≈ 0.57–0.60 y viscosidad ≈ 5,500–9,300 cSt
  al final de corridas de 22–30 días (se acerca al umbral de 10,000 cSt sin
  cruzarlo).

Efectos cuantitativos principales (véase `findings_emulsification_fix_audit24.md`):

| Caso | Evaporado (antes/después) | Biodegradado (antes/después) |
|---|---|---|
| Nohoch 30 d | 85.1 % / 94.1 % | 5.36 % / 0.98 % |
| Abkatun 22 d (siembra superficial) | 81.9 % / 91.5 % | 8.02 % / 1.01 % |
| Abkatun inicio TAMOC (d₀=0.56 mm) | 61.2 % / 74.6 % | 32.0 % / 19.5 % |

El mecanismo: antes de la corrección, la emulsificación rápida encerraba la
mayor parte de la masa de crudo en una emulsión con mucha agua, lo que suprimía
la evaporación y dejaba más crudo disponible para biodegradación. La corrección
restablece una evaporación físicamente realista como el proceso de destino
dominante para el crudo MAYA liberado en superficie en la Bahía de Campeche.

### 8.5 Prueba unitaria

`test_emul_constants.py` verifica (a nivel de código fuente, y de clase cuando
OpenDrift es importable) que las seis constantes valen 1.0 / 1000.0 y que la
fórmula de la tasa no cambió.

`test_emulsification_viscosity.py` (4/4 PASS con los valores anteriores
0.3 / 3.0 / 0.0142). **No se actualizó tras la reversión del 2026-10-29:** las
verificaciones 1, 3 y 4 comparan contra el escalamiento anterior, por lo que se
espera que fallen con los valores 1.0 hasta que la prueba fije explícitamente
los valores anteriores:
- el crudo frío (5 °C, viscosidad alta → estabilidad > 1.0) emulsifica más
  rápido que el cálido (30 °C, viscosidad menor → estabilidad < 1.0)
- con `EMUL_STABILITY_MIN=MAX=1.0`, frío==cálido (se recupera la independencia
  de la viscosidad)
- la estabilidad por defecto siempre es < estabilidad unitaria (el escalamiento
  suprime la tasa para MAYA fresco)
- `EMUL_RATE_SCALE=1.0` da la tasa original de NOAA (mucho más rápida que la
  por defecto)

---

## 8A. Modo de dispersión y tasa de Roy (`processes:dispersion_mode`)

El sumidero de masa de Roy (`processes:dispersion`) y el arrastre de Li et al.
(2017) (`drift:vertical_mixing`) actúan sobre la misma fracción de oleaje roto,
por lo que el crudo puede contarse dos veces. `processes:dispersion_mode`
selecciona cuál actúa:

| Modo | Sumidero de Roy | Arrastre de Li |
|---|---|---|
| `legacy` (por defecto) | activo si `processes:dispersion` es True | activo |
| `transport` | apagado | activo: las gotas entran a la columna de agua y pueden reaflorar |
| `budget` | activo si `processes:dispersion` es True | apagado (estilo PyGNOME/ADIOS) |

El valor por defecto deja sin cambio las configuraciones existentes.

Cambios en el propio sumidero de Roy (2026-10):
- **Solo superficie.** Actúa sobre elementos con `z == 0`; el crudo ya arrastrado
  (`z < 0`) deja de ser removido.
- **Tasa de PyGNOME.** La constante usaba `sqrt(viscosidad / densidad)`,
  dividiendo una viscosidad cinemática entre una densidad. Ahora es
  `sqrt(viscosidad cinemática)` como en `natural_dispersion.disperse_oil_DS` de
  PyGNOME: para 9000 cSt el factor pasa de 0.80 a 9e-4. La fracción es
  `C_Roy · c_disp · V_entrain · dt / (emulsion_density · thickness)`, acotada a 1
  (función `roy_dispersed_fraction()`); el espesor es el `oil_film_thickness`
  del elemento. El crudo viscoso se dispersa mucho menos.
- **No reproducible.** Las corridas con `processes:dispersion` True no
  reproducen la salida anterior a este cambio (corrida de humo, 23 elementos:
  masa dispersada 112,971 → 10,434 kg).

Pruebas: `test_dispersion.py` (20 verificaciones: forma de la tasa, escalamiento,
tope, respaldo de espesor, los elementos sumergidos conservan su masa, seis
combinaciones de modo × `processes:dispersion`).

---

## 8B. Espesor de mancha en malla (`processes:slick_thickness_model = field`)

Por defecto cada elemento lleva su propio espesor (`spillet_thickness`,
esparcimiento de Nihoul si `processes:spreading` es True). En modo `field` la
evaporación lee el espesor de una mancha en malla (`slick_field.py`, portado de
PARCELS_GPU) que es advectada por la deriva superficial y esparcida por la
ecuación del medio poroso, `dh/dt + div(h u) = D0 laplacian(h²)`, con `D0`
calibrado al régimen II gravedad-viscoso de Fay. En cada paso:

1. el campo se escala por B/A, el volumen de emulsión de los elementos que
   estaban en la superficie en la llamada anterior, ahora y entonces (el crudo
   arrastrado y el varado salen como una pérdida uniforme);
2. se advecta con la deriva superficial (corriente + arrastre del viento,
   remuestreada cada `slick_field:drift_update_interval`) y se esparce;
3. el volumen de los elementos que llegan a la superficie (liberación,
   reaflora) se agrega como lentes de Barenblatt (`slick_field:lens_age`);
4. cada elemento superficial dentro de la malla lee el espesor de su celda
   (`spillet_thickness`, y `spillet_area` = volumen de crudo / espesor).

El campo guarda volumen de emulsión por área. `slick_field:evaporation_thickness`
`oil` (por defecto) da a la evaporación la capa de crudo, `h · V_oil /
V_emulsion`; `emulsion` conserva el comportamiento de la primera versión (E1,
1000 m, día 30: evaporado 27.4 % con `emulsion`, 28.5 % con `oil`).

Restricciones y límites:
- Requiere `processes:spreading` False y `processes:julia_weathering` False (el
  modelo lanza `ValueError` en caso contrario).
- La caja de la malla es por defecto el Golfo de México (`lon` −98..−87, `lat`
  15..24.5) con `dx` = 500 m; fije los bordes para otra región. El crudo fuera de
  la caja conserva su propio espesor.
- No implementado: varado en costa en el campo (el volumen que llega a celdas de
  tierra permanece en él), deriva de Stokes en la velocidad de advección, uso del
  campo por la dispersión, la emulsificación y el arrastre (conservan
  `oil_film_thickness`).
- Acoplamiento en un solo sentido: el campo no modifica las partículas, salvo el
  espesor que lee la evaporación.

Pruebas: `test_slick_field.py` (11 verificaciones: conservación de volumen, radio
del frente a 0.7–1.3 % de Fay, traslación, espesor pico de la lente),
`test_slick_coupling.py` (9), `test_slick_lens_thickness.py` (10).

---

## 9. Pruebas y verificación

| Prueba | Qué verifica | Resultado |
|---|---|---|
| `test_cicoil_e2e.py` | Batería completa de 30 pruebas (T1–T8): siembra, intemperización, siembra de pluma TAMOC, balance de masa/componentes, biodegradación (T7), fotooxidación + seguimiento de OP + cota de regresión (T8) | 30/30 PASS, conservación de masa 100.000 %, error de balance por componente 0.000 % |
| `test_photoox_vertmix.py` | `vertical_mixing` activado/apagado contra la referencia de `surface_fate_coupled.py` | `vertical_mixing=True`: 3.50 % contra objetivo 3.47 % (dentro de 1 %); `False`: sobrestima ~11 % |
| `test_biodegradation.py` | Prueba de cordura de Q10 (liberación superficial, 10 d, MAYA, resto de la intemperización apagada) | Tras recalibrar: 0.273 %/0.810 %/2.357 % biodegradado a 10/20/30 °C (contra 7.110 %/16.595 %/29.301 % antes de recalibrar, ~30× menor como se espera); conservación de masa exacta |
| `test_biodegradation_water.py` | Decaimiento/desactivación del reservorio disuelto (20 partículas, `is_dissolved=1`, 60 d, SST 10/20/30 °C) | 10 °C: 68.526 % biodeg, 0/20 desactivadas; 20 °C: 96.882 %, 0/20; 30 °C: 99.043 %, 20/20 desactivadas; conservación de masa ~1e-15 |
| `test_droplet_size_biodeg.py` | Biodegradación resuelta por tamaño de gota (auditoría #2; 90 partículas, 4 grupos de diámetro, 10 d, SST=20 °C) | `frac_biodeg(0.1mm)=7.110% > frac_biodeg(1mm)=0.810% > frac_biodeg(1cm)=0.082%`; `diameter=0` == `diameter=1mm` (dif. rel. 4.7e-8); error relativo de conservación de masa 1.22e-16 |
| `test_op_tracking.py` | Seguimiento de fotoproductos oxigenados (auditoría #4; 100 partículas, 8 d, SST=29 °C, UV sintético) | conservación de masa incl. reservorios OP, error rel. 1.22e-16; `mass_op_dissolved + mass_op_degraded == mass_photooxidized` (OP_DISSOLVED_FRACTION=0.5, exacto); `mass_op_degraded > 0` (decaimiento K_OP_REMOVAL activo) |
| `test_emulsification_viscosity.py` | Escalamiento viscosidad-estabilidad: el crudo frío emulsifica más rápido que el cálido; estabilidad unitaria recupera independencia de la viscosidad; la estabilidad por defecto suprime la tasa frente a la unitaria; EMUL_RATE_SCALE=1.0 recupera la tasa original de NOAA | 4/4 PASS con los valores anteriores; requiere actualización (§8.5) |

| `test_horizontal_diffusivity.py` | Paseo aleatorio con forzamiento constante, 6 h, 2000 partículas: D=0 no dispersa; D=10 m² s⁻¹ vía `drift:horizontal_diffusivity` o vía el respaldo del ambiente se dispersa como √(2Dt) (659 m contra 657 m); D=40 se dispersa 2.00× D=10 | 4/4 PASS |
| `test_dispersion.py` | Forma de la tasa de Roy, escalamiento por espesor/densidad, tope en 1, elementos sumergidos intactos, seis combinaciones dispersion_mode × `processes:dispersion` | 20 verificaciones |
| `test_slick_field.py` / `test_slick_coupling.py` / `test_slick_lens_thickness.py` | Mancha en malla: frente de Fay, conservación, acoplamiento con los elementos, espesor de crudo frente a emulsión | 11 / 9 / 10 verificaciones |
| `test_evaporation_components.py` | La suma de `comp_*_evaporated` es igual a `mass_evaporated`; remanente + evaporado constante; el residuo no se evapora | véase el archivo |
| `test_local_sea_water_density.py` | T y S del elemento usadas para la densidad del agua de mar en la intemperización | véase el archivo |
| `test_ift_units.py` | IFT en N/m (por defecto 0.03) | véase el archivo |
| `test_emul_constants.py` | Las constantes de emulsificación valen 1.0 | véase el archivo |
| `test_croco_native_reader.py` | `reader_croco_native` en un archivo CROCO sintético: `mask_u`/`mask_v` derivadas, `u`/`v` desescalonadas en un nodo rho iguales al valor analítico, `destagger=False` muestra el desfase de media celda original | 7/7 PASS |
| `test_deactivate_masks.py` | Las burbujas de gas liberadas y los elementos disueltos agotados se desactivan con una máscara booleana (OpenDrift 1.14.9 la exige); ambas corridas terminan | requiere la rama `fix/gas-bubble-mask` |
| `test_oil_density_ids.py` | Siembra escalonada de pluma (gas, crudo, disuelto): ningún elemento de crudo queda activo con masa NaN (17 antes de la corrección, 0 después) | requiere la rama `fix/gas-bubble-mask` |
| `test_vertical_profile.py` | `vertical_profile.sample_profile` con lectores sintéticos (columnas crudas, NaN bajo el fondo se conserva, se agrega el nivel superficial); la prueba del escritor TAMOC se omite si TAMOC no está instalado | requiere la rama `feat/vertical-profile` |

### Advertencias conocidas
- **Emulsificación**: se requiere `max_water=0.9` en `seed_elements()` para MAYA
  (ADIOS AD02254 no tiene datos de emulsificación → de otro modo
  `emulsification_cicese()` termina de inmediato). Desde el 2026-10-29 la tasa es
  la original de NOAA (`EMUL_*` = 1.0); la calibración anterior de bolas de
  alquitrán a 17.2 días se retiró (§8).
- **ANFC frente a GLORYS12**: el forzamiento oceánico ANFC da resultados no
  reproducibles entre fechas de descarga (80 km contra 151 km de transporte para
  la misma corrida nominal). Use GLORYS12 (`--ocean glorys`) para atribución
  reproducible en hindcast.
- **Signo de la longitud**: la Bahía de Campeche es °W: use `lon=-92.x`, no
  `+92.x` (una longitud positiva coloca las partículas en el océano Índico sin
  aviso).
- **Masa por elemento**: CICOIL no impone una masa fija por elemento; es la
  masa de crudo sembrada dividida entre el número de elementos. El valor
  185.6 kg heredado de v1 no lo aplica el modelo (el script TAMOC de Nohoch 2023
  lo usó solo para escalar gráficas, frente a 0.31 kg reales por elemento).
  Escale las salidas con `M_spill / (N_p × masa por elemento)` usando la masa
  realmente sembrada (suma de `mass_oil` en el primer paso de salida).

---

## 10. Manifiesto de archivos (paquete de despliegue en HPC)

**Modelo principal:**
- `ciceseoil.py` (142 KB) — clase `OpenCiceseOil`, 14 correcciones de API, fotooxidación + biodegradación + emulsificación
- `slick_field.py` — campo de espesor de mancha en malla (§8B)
- `fluid_properties.py` (19 KB) — marco de destilación de 17 componentes (Pv de Cox + Mackay)
- `cicoil_estimations.py` (13 KB) — soluciones independientes de estimación de la base ADIOS

**Acoplamiento con TAMOC:**
- `tamoc_plume.py`, `run_tamoc.py`, `tamoc_chemical_properties.py`, `vertical_profile.py` (rama `feat/vertical-profile`, §10C)

**Lectores NEMO/CROCO** (`readers/`):
- `reader_NEMO_native_v3.py`, `reader_nemo_combined.py`, `reader_nemo_modified.py`, `reader_nemo_optimized.py`
- `reader_croco_native.py` — salida cruda de CROCO en coordenada sigma, desescalonada (subclase de `reader_ROMS_native` de OpenDrift)

**Exportación:**
- `export/io_stat_nc.py` — salida NetCDF estadística en malla

**Datos** (`data/`, 5 CSV, ~200 KB):
- `pseudo_chemdata.csv`, `ChemData.csv`, `Aij.csv`, `Bij.csv`, `gas_composition.csv`

**Pruebas:**
- `test_cicoil_e2e.py` (batería de 30 pruebas, T1–T8), `test_biodegradation.py`,
  `test_biodegradation_water.py`, `test_droplet_size_biodeg.py`,
  `test_op_tracking.py`, `test_emulsification_viscosity.py`

**Auxiliares compartidos:**
- `cicoil_common.py` — `load_ciceseoil()` (código de importlib) y
  `default_cicese_config()` (configuración estándar de procesos activados/apagados)

**Utilidades:**
- `download_data.py` (auxiliares de descarga de GLORYS12/ERA5)
- `preprocess_era5_uv.py`, `preprocess_real_forcing.py`
- `plotting.py` (10+ rutinas de graficación)

**Despliegue:**
- `CICOIL_v2_MANUAL.md` (este archivo, en inglés), `stage_for_hpc.sh`, `deploy_on_hpc_conda.sh`,
  `slurm_template_cicoilv2.sh`, `DEPLOYMENT_MANIFEST.txt`

Total: ~500 MB el paquete, ~150 MB instalado en site-packages.

---

## 10A. Cambios desde 28f8cdb (2026-10-08), rama `physics/e1-weathering-fixes`

| Cambio | Efecto en los resultados |
|---|---|
| Kelvin/Celsius en `diffusion_coeff_in_air` (`+ 273.15`) | La evaporación era 25–30 % demasiado baja en agua cálida (Abkatun adelante, 47 d: 0.1583 contra 0.2059 de la masa liberada) |
| T y S del elemento para la densidad del agua de mar (`local_sea_water_density()`) | La densidad de la emulsión era unos 3.2 kg m⁻³ demasiado alta (30.2 °C, S 36.6: 1022.86 contra 1026.95 kg m⁻³) |
| Dispersión de Roy solo en elementos superficiales; tasa `sqrt(nu)` (§8A) | Mucha menos dispersión de crudo viscoso; las corridas con `processes:dispersion` True cambian |
| `processes:dispersion_mode` (§8A) | Por defecto `legacy`: sin cambio |
| Campo de mancha en malla y `processes:slick_thickness_model` (§8B) | Por defecto `legacy`: sin cambio |
| `comp_*_evaporated` guarda la masa evaporada acumulada (antes la masa remanente) | Solo las salidas por componente; los balances no cambian (variables de solo escritura) |
| IFT en N/m, por defecto 0.03 (antes 30) | La tasa de arrastre de Li et al. (2017) y el espectro de gotas usan N/m |
| Constantes `EMUL_*` de vuelta a 1.0 (§8) | Las corridas sin ajuste en tiempo de ejecución emulsifican más rápido |

---

## 10B. Migración desde CICOIL v1

Para quienes usan CICOIL v1 (`Cicoil_Version_22-06-2023`, un fork de OpenDrift
1.9.0, TAMOC 3.4.0). v1 está congelado y no se corregirá; todo el desarrollo
está en CICOILv2 (OpenDrift 1.14.9, TAMOC 4.1.2).

### Qué cambia para el usuario

| Tema | v1 | CICOILv2 |
|---|---|---|
| Base | fork de OpenDrift 1.9.0 | OpenDrift 1.14.9 sin modificar; los archivos de CICOIL se copian en `opendrift/models/openoil/` |
| TAMOC | 3.4.0 | 4.1.2 |
| Partículas perdidas | guardadas como `inf` / 1e37 | guardadas como NaN: use `np.nanmean` / `np.nansum` o filtre por `status` |
| Números aleatorios | — | pase `seed=` al constructor del modelo; OpenDrift controla su generador |
| Claves de configuración | el lector de v1 rechaza claves desconocidas | claves de OpenDrift 1.14.9, p. ej. `drift:profiles_depth` (plural); copie una configuración antigua a un script nuevo y deje que OpenDrift señale las claves inválidas |
| Lector de ERA5 | falla con las coordenadas escalares `number` / `expver` | las lee (v1 necesitaba un ERA5 preprocesado) |
| Salida de CROCO | requiere conversión a niveles z en malla regular | `reader_croco_native` lee la salida cruda en coordenada sigma |
| Fotooxidación, biodegradación | no disponibles | opcionales, apagadas por defecto (§6, §7) |

### Por qué difieren los resultados de v1 y CICOILv2

Tres errores encontrados al correr v1 contra OpenOil, PyGNOME y CICOILv2 en
Nohoch 2023 cambian los resultados de v1 (informe `V1_BUGS` en
`CICOIL_COMPARISON`):

1. **Kelvin/Celsius** en `diffusion_coeff_in_air`: la evaporación queda demasiado
   baja (v2: 25–30 % en agua cálida).
2. **Velocidad terminal de gotas**: la temperatura en K entra en la densidad del
   agua de mar y el diámetro se usa donde se necesita un radio, de modo que las
   gotas arrastradas regresan a la superficie casi de inmediato y v1 muestra
   casi nada de crudo sumergido.
3. **Constante de dispersión de Roy** usa `sqrt(nu/rho)` en vez de `sqrt(nu)`
   (la misma expresión está en OpenDrift 1.14.9 original): E1 dispersa 92 % con
   la forma programada contra 3.3 % con la forma de PyGNOME.

CICOILv2 corrige los tres. Por lo tanto, no se espera que una corrida de v1 y
una de v2 con la misma configuración coincidan, y la diferencia no es ruido.

### Instalación compartida en Ixachi (CICESE)

Todo en `/LUSTRE/adomingu/CICOIL_SHARED` (solo lectura para los usuarios;
copie a su espacio lo que necesite):

- `env/cicoilv2`: entorno conda con OpenDrift 1.14.9, CICOILv2 (commit
  `b032140`) y TAMOC 4.1.2 (copia privada parchada en `src/tamoc_4.1.2`; con
  numpy 2.x hace falta `f2py --backend meson`). Los metadatos del paquete TAMOC
  dicen 4.0.0; `tamoc.__version__` dice 4.1.2.
- `examples/abkatun2026/inputs` (conjunto pequeño GLORYS + ERA5) y
  `examples/nohoch2023`: forzamiento preprocesado.
- `ejemplo/abkatun2026.py` con `run_abkatun2026.slurm` (`MODE=2d|3d`,
  `FORCING=coupled|glorys`, `DROPLETS=tamoc|default|lognormal:D50:SIGMA`);
  `ejemplo/run_nohoch2023.slurm`.
- `readers_check/` y `hpc/10_check_readers.slurm`: pruebas de lectores para
  CROCO (nativo), WW3 y WRF (mediante archivos CF creados con `prep_forcing.py`;
  OpenDrift no tiene lector nativo de WW3/WRF) y NEMO (se omite hasta dar un
  archivo NEMO).
- Nada de cómputo en el nodo de acceso: todo con `sbatch` (partición
  `oc_cpus`; `oc_gpu` para GPU).
- Entorno propio: `deploy_on_hpc_conda.sh` (§2.2).

### Trampas

- **Longitud**: el Golfo es oeste, `lon=-92.x`; `+92` coloca las partículas en
  el océano Índico sin aviso.
- **`max_water`**: sin `max_water=0.9` en `seed_elements()`, MAYA no emulsifica.
- **Reproducibilidad**: las corridas hechas con `processes:dispersion` True
  antes de octubre de 2026 no se reproducen (§8A).
- Guarde la semilla y el commit con cada resultado.

### Estado de validación

- S0 (Nohoch 2023, día 30, `EMUL_*` = 1.0): modo `budget` 2.1 % disperso /
  26.9 % evaporado / 71.0 % remanente; `legacy` 1.7 / 23.8 / 74.6; PyGNOME
  4.4 / 24.4 / 70.6. `budget` queda 2.3 puntos por debajo de PyGNOME en
  dispersión; la causa rastreada es el momento en que empieza la emulsificación
  (CICOIL fija 2 h; PyGNOME usa la fracción de Bullwinkle del registro del
  crudo), no el espesor de la mancha.
- Pendientes: criterio basado en SAR; `dispersion_mode` por defecto para E1;
  escala de masa por elemento de PyGNOME; campos subsuperficiales de
  `reader_croco_native`.

---

## 10C. Cambios posteriores a b032140 (ramas aún sin fusionar)

| Rama | Cambio | Efecto en los resultados |
|---|---|---|
| `fix/gas-bubble-mask` | `handle_gas_particles` y `handle_subsea_dissolution` pasan una máscara booleana a `deactivate_elements` (OpenDrift 1.14.9 da un error de broadcast con índices enteros) | Las corridas con burbujas de gas (`is_bubble=1`) o elementos disueltos agotados se detenían en el primer paso; las demás no cambian |
| `fix/gas-bubble-mask` | `oil_density` indexa `cicese_mass_balance` por ID de elemento, como `component_molfrac` (usaba posiciones del arreglo; v1 usaba `indices - 1`) | Idéntico mientras todos los elementos se siembran juntos y ninguno se elimina. Con siembra escalonada (plumas TAMOC) o elementos eliminados, elementos de crudo leían la fila de masa de otro elemento, quedaban con status 0 y densidad y masa NaN y salían del balance (Abkatun 2026 3D: 119 de 427 elementos, cierre 67 %). Las corridas en que se eliminaron elementos a mitad (varamiento, salida del dominio) no se han auditado |
| `feat/vertical-profile` | Nuevo `vertical_profile.py`: `create_profile()` construye el perfil ambiental de TAMOC (T, S, u, v con la profundidad) a partir de los lectores de OpenDrift; `sample_profile()` devuelve las columnas | Necesario para la ruta TAMOC siguiente; no afecta otras corridas |

### Ruta TAMOC (como en v1)

1. `get_live_composition` para crudo y gas;
2. `vertical_profile.create_profile(nc_file, readers, lon, lat, time, ...)`
   escribe el perfil ambiental; lee la columna cruda de la celda más cercana con
   `reader.get_variables`, no con `get_variables_interpolated`, que rellena los
   NaN (fondo, tierra) con vecinos y devuelve bloques verticales viejos si se
   llama nivel por nivel;
3. `run_tamoc` ejecuta el modelo de pluma curva;
4. `Plume(tracked=False)`: TAMOC 4.1.2 no escribe archivos de trayectoria por
   partícula, por lo que `tracked=True` falla;
5. `seed_plume_elements(seed_time=[inicio, fin])`.

El forzamiento debe tener niveles verticales: los archivos GLORYS/ANFC de
Abkatun 2026 tienen un nivel (0.49 m); el ejemplo 3D usa los archivos acoplados
CROCO-WRF-WW3 (10 niveles). `create_profile(..., z_max=profundidad+20)` extiende
la columna por debajo de la profundidad de liberación.

---

## 11. Referencias

- Adcroft et al. (2010) — ley de escalamiento por temperatura Q10=3, T_REF=20 °C
  para biodegradación; tasa de referencia de pluma disuelta del Golfo de México
- Bacosa et al. (2015), Marine Pollution Bulletin — rangos de vida media por
  clase de anillo para biodegradación/fotooxidación
- Prince et al. (2003, 2017), ACS ES&T "The Rate of Crude Oil Biodegradation in
  the Sea" — contraste de biodegradación entre mancha superficial y crudo
  disperso (base de `BIODEG_SURFACE_FACTOR`)
- Delvigne & Sweeney (1988), Oil & Chem. Pollut. 4:281 — sumergimiento por oleaje
  (histórico; sustituido por la fracción superficial basada en `vertical_mixing`)
- Ward et al. (2018) — fotooxidación alifática mediada por radicales (discutida,
  no implementada: solo aromáticos es una cota inferior)
- Documentación de OpenDrift 1.14.9
- Informes técnicos de CICOIL 2023 (Kotzakoulakis, CICESE/SINTEF Ocean)

**Bitácoras de sesión:**
- `session_2026-06-05_CICOIL_analysis.{html,pdf}` — auditoría completa del código
- `session_2026-06-06_CICOIL_install_test.{html,pdf}` — instalación y pruebas
- `session_2026-06-13_era5_uv_recalibration.{html,tex}` — recalibración del UV de fotooxidación
- `session_2026-06-14_CICOILv2_Biodegradation.{tex,pdf}` (12 pp) — módulo de biodegradación
- `BIODEGRADATION.md` — diagnóstico completo de biodegradación + tablas de tasas

---

**Versión:** CICOILv2 (commit `b032140` más las ramas de la §10C; etiqueta de versión `v2.1.0` pendiente)
**Última actualización:** 2026-11-17
