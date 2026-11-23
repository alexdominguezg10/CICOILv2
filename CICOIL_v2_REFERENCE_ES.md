# Manual de referencia de CICOILv2

> Traducción al español de `CICOIL_v2_REFERENCE.md` (la versión en inglés es la referencia; si difieren, manda el inglés). Los nombres de funciones, claves y variables no se traducen.

**Referencia completa de la funcionalidad de `cicoilv2_deployment/`**

Este es el documento de referencia complementario de
[`CICOIL_v2_MANUAL_ES.md`](CICOIL_v2_MANUAL_ES.md) (resumen, despliegue,
calibración de fotooxidación/biodegradación). Mientras el Manual explica *cómo
desplegar y correr* CICOIL v2 y la física de las adiciones de intemperización de
CICESE, esta Referencia cataloga *cada módulo, clase, método, función, variable
de elemento y clave de configuración* disponible en el paquete, para consulta de
la API.

**Fuente:** `OILSPILL/CICOIL_dev/cicoilv2_deployment/`

---

## 1. Estructura del paquete

```
cicoilv2_deployment/
├── ciceseoil.py              # CiceseOil (tipo de elemento) + OpenCiceseOil (clase del modelo)
├── fluid_properties.py       # FluidProps: marco de destilación de 17 componentes
├── cicoil_estimations.py      # Funciones independientes de estimación de la base ADIOS
├── slick_field.py             # campo de espesor de mancha en malla (processes:slick_thickness_model=field), véase §3.3
├── tamoc_plume.py             # Clase Plume — salida BPM de TAMOC → siembra en OpenDrift
├── run_tamoc.py                # run_tamoc() — envoltura del BPM de TAMOC
├── tamoc_chemical_properties.py
├── readers/
│   ├── reader_NEMO_native_v3.py    # activo — nombres de variable nativos de CROCO
│   ├── reader_nemo_optimized.py    # activo — nombres genéricos NEMO-CMEMS
│   ├── reader_nemo_combined.py     # alternativo — basado en xarray, véase §8.3
│   ├── reader_nemo_modified.py     # alternativo — basado en netCDF4, véase §8.3
│   └── reader_croco_native.py      # activo — salida cruda de CROCO en coordenada sigma, véase §8.4
├── export/
│   └── io_stat_nc.py           # exportación NetCDF estadística en malla
├── plotting.py                 # auxiliares de graficación para publicación (8 funciones)
├── download_data.py            # descarga y preprocesamiento de GLORYS12/ANFC/ERA5
├── preprocess_era5_uv.py        # ERA5 UV → NetCDF conforme a CF (forzamiento de fotooxidación)
├── preprocess_real_forcing.py
├── cicoil_common.py             # auxiliares compartidos (load_ciceseoil, default_cicese_config)
├── test_cicoil_e2e.py           # batería completa de 30 pruebas (T1–T8)
├── test_biodegradation.py        # prueba de cordura de Q10
├── test_biodegradation_water.py  # prueba de la vía de fase disuelta
├── CICOIL_v2_MANUAL.md           # despliegue + física de los módulos de intemperización
└── CICOIL_v2_REFERENCE.md         # este archivo
```

---

## 2. Variables de elemento (`CiceseOil`)

`CiceseOil` extiende el tipo de elemento `Oil` de OpenDrift. Total: **97
variables** (22 heredadas de `Oil`/`OpenOil` original + 75 específicas de CICOIL
abajo).

### 2.1 Variables de estado centrales de CICOIL

| Variable | Unidades | Valor por defecto | Descripción |
|---|---|---|---|
| `spillet_area` | m² | 0.114 | Área superficial de una cara del spillet |
| `oil_molar_mass` | kg/mol | 0.2 | |
| `spillet_thickness` | m | 0.01 | |
| `oil_density` | kg/m³ | 880 | |
| `DOR` | — | 0 | Razón dispersante-crudo (Dispersant-to-Oil Ratio) |
| `IFT` | N/m | 0.03 | Tensión interfacial crudo–agua (N/m; valía 30.0 en mN/m antes del 2026-10-26, lo que daba a la tasa de Li et al. 2017 y al espectro de gotas un valor 1000× demasiado grande) |
| `max_water` | — | 0 | Contenido máximo de agua (tope de emulsificación) |
| `mass_dissolved_surface` | kg | 0 | |
| `mass_dissolved_subsea` | kg | 0 | Se usa como total de referencia del reservorio disuelto (§7.3 del Manual) |
| `is_bubble` | — | 0 | Bandera: partícula de burbuja de gas |
| `is_dissolved` | — | 0 | Bandera: partícula disuelta submarina |
| `mass_gas` | kg | 0 | |
| `mass_biodegraded_from_oil` | kg | 0 | NUEVO — la escribe `biodegradation_cicese()` |
| `mass_biodegraded_from_water` | kg | 0 | NUEVO — la escribe `biodegradation_cicese()` (vía disuelta) |
| `fraction_dissolved_surface` | % | 0 | |
| `fraction_biodegraded_from_oil` | % | 0 | NUEVO |
| `fraction_biodegraded_from_water` | % | 0 | NUEVO |
| `mass_photooxidized` | kg | 0 | NUEVO — la escribe `photooxidation_cicese()` |
| `fraction_photooxidized` | % | 0 | NUEVO |
| `mass_op_dissolved` | kg | 0 | NUEVO (auditoría #4) — reservorio disuelto de fotoproductos oxigenados, alimentado por `photooxidation_cicese()` |
| `fraction_op_dissolved` | % | 0 | NUEVO (auditoría #4) |
| `mass_op_degraded` | kg | 0 | NUEVO (auditoría #4) — producto del decaimiento de `mass_op_dissolved` vía `K_OP_REMOVAL` |
| `fraction_op_degraded` | % | 0 | NUEVO (auditoría #4) |

### 2.2 Arreglos de pseudo-componentes (17 componentes × 3 = 51 variables)

Para cada uno de los 17 pseudo-componentes siguientes existen tres variables por
partícula: `<comp>`, `<comp>_evaporated`, `<comp>_dissolved` (todas en kg, valor
por defecto 0).

```
comp_ali_G1 .. comp_ali_G8     (8 grupos alifáticos)
comp_arom_G1 .. comp_arom_G8   (8 grupos aromáticos)
comp_residue                   (1 residuo no volátil)
```

G1 = el más ligero (p. ej. benceno/alcanos ligeros) → G8 = el más pesado
(pireno+/ceras pesadas). Las leen y escriben `evaporation_cicese()`,
`biodegradation_cicese()` y `photooxidation_cicese()` (solo aromáticos), y se
siguen en `self.cicese_mass_balance['mass_components']`.

### 2.3 Variables heredadas (de `OpenOil`/`Oil` original, selección)
`mass_oil`, `mass_evaporated`, `mass_dispersed`, `mass_biodegraded`,
`density`, `viscosity`, `water_fraction`, `age_seconded`, `oil_film_thickness`,
`diameter`, `wind_drift_factor`, `bulltime`, más posición/estado/ID — véase la
documentación de `OpenOil` de OpenDrift para el conjunto completo de 22
variables heredadas.

---

## 3. Claves de configuración

### 3.1 Claves específicas de CICOIL/CICESE y agregadas en OD 1.14.9

| Clave | Tipo | Valor por defecto | Descripción |
|---|---|---|---|
| `seed:m3_per_hour` | float | 1 | Volumen de crudo liberado por hora (o total si es instantáneo) |
| `seed:droplet_diameter_min_subsea` | float | 0.0005 m | Diámetro mínimo de gota, liberación submarina |
| `seed:droplet_diameter_max_subsea` | float | 0.005 m | Diámetro máximo de gota, liberación submarina |
| `seed:droplet_size_distribution` | enum | `'uniform'` | `uniform`/`normal`/`lognormal`, liberación submarina |
| `seed:droplet_diameter_mu` | float | 0.001 m | Diámetro medio para normal/lognormal |
| `seed:droplet_diameter_sigma` | float | 0.0005 m | Desviación estándar para normal/lognormal |
| `seed:oil_type` | enum | el primero de `self.oiltypes` | Nombre del crudo en ADIOS |
| `processes:dispersion` | bool | False | Dispersión natural (sumidero de masa de Roy; véase `processes:dispersion_mode`) |
| `processes:dispersion_mode` | enum | `legacy` | `legacy` / `transport` / `budget`: cuál de los dos, el sumidero de masa de Roy o el arrastre de Li et al. (2017), actúa. Véase §3.3 y Manual §8A |
| `processes:slick_thickness_model` | enum | `legacy` | `legacy` / `field`: espesor de mancha que lee la evaporación. `field` requiere `processes:spreading` False y `processes:julia_weathering` False. Véase §3.3 y Manual §8B |
| `processes:chemical_dispersion` | enum | No\_treatment | No\_treatment / Surface / SSDI |
| `processes:evaporation` | bool | True | Evaporación por Pv de Cox |
| `processes:emulsification` | bool | True | Captación de agua basada en Bullwinkle |
| `processes:biodegradation` | bool | **False** | NUEVO — `biodegradation_cicese()`, véase Manual §7 |
| `processes:update_oilfilm_thickness` | bool | False | Recalcular el espesor de la película en cada paso o mantenerlo constante |
| `processes:handle_released_gas` | bool | True | Contabilidad de burbujas de gas (`handle_gas_particles()`) |
| `processes:subsea_dissolution` | bool | True | Alimenta el reservorio de masa disuelta (`handle_subsea_dissolution()`) |
| `processes:spreading` | bool | True | Esparcimiento de Nihoul (1984) |
| `processes:surface_dissolution` | bool | False | |
| `processes:photooxidation` | bool | **False** | NUEVO — `photooxidation_cicese()`, véase Manual §6 |
| `wave_entrainment:` `droplet_size_distribution` | enum | Johansen et al. (2015) | Johansen et al. (2015) / Li et al. (2017) |
| `wave_entrainment:entrainment_rate` | enum | Li et al. (2017) | |
| `biodegradation:method` | enum | `'Adcroft'` | Forzado a `'Adcroft'`; describe la ley Q10 de `biodegradation_cicese()` (no la fórmula global original) |

### 3.2 Valores de deriva fijados por `OpenCiceseOil.__init__`

| Clave | Valor | Notas |
|---|---|---|
| `drift:vertical_advection` | False | |
| `drift:vertical_mixing` | True | **Requerido** para la fracción superficial correcta en fotooxidación (Manual §6.2) |
| `drift:current_uncertainty` | 0.05 | |
| `drift:wind_uncertainty` | 0.5 | |
| `drift:profiles_depth` | 20 | Clave en plural de OD 1.14.9; sustituye al obsoleto `required_profiles_z_range` |

---

### 3.3 Claves del modo de dispersión y del campo de mancha

`processes:dispersion_mode` (enum, por defecto `legacy`). El sumidero de masa de Roy (`processes:dispersion`) y el arrastre de Li et al. (2017) (`drift:vertical_mixing`) actúan sobre la misma fracción de oleaje roto, por lo que `legacy` cuenta dos veces el mismo crudo.

| Modo | Sumidero de Roy | Arrastre de Li |
|---|---|---|
| `legacy` | activo si `processes:dispersion` es True | activo |
| `transport` | apagado, diga lo que diga `processes:dispersion` | activo (gotas transportadas, pueden reaflorar) |
| `budget` | activo si `processes:dispersion` es True | apagado (estilo PyGNOME/ADIOS) |

Claves `slick_field:*` (se leen solo cuando `processes:slick_thickness_model` = `field`):

| Clave | Tipo | Valor por defecto | Descripción |
|---|---|---|---|
| `slick_field:dx` | float | 500 m | Tamaño de celda de la malla de espesor |
| `slick_field:lon_min` / `lon_max` | float | -98 / -87 | Borde oeste / este de la malla (°E); el crudo fuera de la caja no es seguido por el campo |
| `slick_field:lat_min` / `lat_max` | float | 15 / 24.5 | Borde sur / norte de la malla (°N) |
| `slick_field:lens_age` | float | 300 s | Edad de la lente de Barenblatt en la que se deposita el volumen de un elemento recién llegado a la superficie |
| `slick_field:evaporation_thickness` | enum | `oil` | `oil`: la evaporación lee h·V_oil/V_emulsion (volumen de crudo por área, como el esparcimiento por elemento); `emulsion`: el espesor de la celda tal cual (comportamiento del commit 66bce53) |
| `slick_field:drift_update_interval` | float | 3600 s | Intervalo al que se remuestrea en la malla la deriva superficial (corriente + arrastre del viento) |

La caja de malla por defecto es el Golfo de México: fije los cuatro bordes para cualquier otra región.

`slick_field.py` (NumPy, portado de PARCELS_GPU): `Grid`, `fay_D0()`, advección MUSCL/minmod con pasos de Heun, esparcimiento por medio poroso, depósito de lente de Barenblatt, diagnósticos de clase Bonn; auxiliares `emulsion_volume()` y `sample()`.

## 4. `OpenCiceseOil` — referencia de métodos

`class OpenCiceseOil(OpenOil)` en `ciceseoil.py`. Tipo de elemento:
`ElementType = CiceseOil`. `required_variables` incluye los campos estándar de
corriente/viento/oleaje/hielo/temperatura/salinidad del océano más
`surface_downward_uv_radiation` (respaldo `-1.0`, usado por
`_compute_uv_irradiance()`).

### 4.1 Construcción y selección del crudo

| Método | Firma | Descripción |
|---|---|---|
| `__init__` | `(self, weathering_model='noaa', *args, **kwargs)` | `weathering_model` debe ser `'noaa'` o `'cicese'`. Construye `self.oiltypes` a partir de la base ADIOS (primero los crudos GENERIC, filtrando `duplicate_oils`), registra todas las claves de configuración (§3), fija los valores de deriva. |
| `set_oiltype` | `(self, oiltype)` | Fija `self.oiltype` (un objeto fluido `Oil`/`CiceseOil` de `adios_db`) y, si `weathering_model='cicese'`, construye `self.fluid_properties = FluidProps(...)` (§5) para el marco de 17 componentes. |
| `get_max_water_content` | `(self)` | Devuelve la fracción máxima de agua (emulsificación) del crudo, aplicando `_cicoil_max_water_fraction_override` para los dos crudos con vacíos conocidos en ADIOS (`MARINE GAS OIL 500 ppm S 2017`, `FENJA (PIL) 2015`). |
| `run` | `(self, record_components=False, record_environment=True, *args, **kwargs)` | Envoltura ligera de `OpenOil.run()`. Por defecto `export_variables` excluye los 51 arreglos por componente (`comp_*`, `comp_*_evaporated`, `comp_*_dissolved`), que de otro modo dominan el tamaño de la salida en ensambles grandes. `record_components=True` los conserva. `record_environment=False` además omite los 20 campos de ambiente de `required_variables` (corriente, viento, oleaje, SST, salinidad, etc.) muestreados por elemento: útil para reducir aún más la salida de ensambles grandes, porque se pueden recuperar de los datos de forzamiento. Ambas banderas se ignoran si se pasa un argumento `export_variables` explícito. OpenDrift comprime la salida NetCDF con zlib (complevel 6) al cerrar el archivo, independientemente de estas banderas. |
| `prepare_run` | `(self)` | Preparación previa a la corrida: inicializa los arreglos de seguimiento por componente de `cicese_mass_balance`, llena `self.max_water_fraction` desde el JSON de ADIOS (modo NOAA). |

### 4.2 Siembra

| Método | Firma | Descripción |
|---|---|---|
| `seed_elements` | `(self, *args, **kwargs)` | Punto de entrada principal de siembra. Acepta `lon`, `lat`, `time`, `number`, `z`, `m3_per_hour`, `oil_type` (preferido) u `oiltype` (obsoleto), `max_water`, `is_bubble`, `is_dissolved`, `diameter`. Calcula `mass_oil` por elemento como `m3_per_hour × duración / number × densidad`; para `z<0` submarino sin `diameter` explícito, extrae diámetros de gota uniformemente entre `seed:droplet_diameter_min/max_subsea`. Fija `self.Density`/`self.KinematicViscosity` (modo NOAA) o `density`/`oil_density`/`viscosity` (modo cicese) a partir de `oiltype.density_at_temp(285)` / `kvis_at_temp(285)`. `is_bubble=1` → la masa se dirige a `mass_gas`; `is_dissolved=1` → la masa se dirige a `mass_dissolved_subsea`. |
| `seed_with_dispersant` | `(self, lon, lat, seed_time, appl_method, DOR, dispersant='C9500', time_frac_surf, time_frac_ssdi, effic_surf, effic_ssdi, wind_limit, delay, stock_kg, time_step, appl_rate, number)` | Siembra elementos que representan una mancha dispersada químicamente: calcula la IFT/max-water tratadas con `treated_IFT_and_WC()`, muestrea el DOR efectivo con `dispersant_efficiency_sampled()` y verifica las ventanas de clima con `weather_permitting_time()`. |
| `seed_plume_elements` | `(self, lon, lat, plume, seed_time, *args, **kwargs)` | Siembra partículas desde un objeto `Plume` de TAMOC (§7.1), es decir una descarga submarina. Distribuye elementos de gotas/burbujas de gas/fase disuelta según la salida de seguimiento de partículas de la pluma (`Plume.get_particle_properties()`), fija las banderas `is_bubble`/`is_dissolved` y finaliza el ambiente (con un respaldo de temperatura de 28 °C si aún no está disponible) antes de pasar el modelo de `Mode.Config` a `Mode.Ready` mediante `__set_seed_config__()`. |
| `weather_permitting_time` | `(self, lon, lat, treatment_time, wind_limit, time_step)` | Calcula la ventana de tiempo acumulada durante la cual la rapidez del viento está por debajo de `wind_limit` en la posición/tiempo dados; se usa para condicionar la aplicación de dispersante. |

### 4.3 Procesos de intemperización (despacho y física)

| Método | Firma | Descripción |
|---|---|---|
| `oil_weathering` | `(self)` | Despachador de nivel superior; llama a `oil_weathering_cicese()` (o a la ruta NOAA original) según `oil_weathering_model`. |
| `oil_weathering_cicese` | `(self)` | Tubería de intemperización por paso de tiempo (esquema CICOIL adaptado de NOAA PyGNOME). Llama, en orden: `evaporation_cicese()`, `emulsification_cicese()`, `spreading_cicese()`, `dispersion_cicese()`, `handle_subsea_dissolution()`, `handle_gas_particles()`, y luego, condicionalmente, `biodegradation_cicese()` (si `processes:biodegradation`) y `photooxidation_cicese()` (si `processes:photooxidation`). |
| `evaporation_cicese` | `(self, atm_pressure=101325.)` | Evaporación por presión de vapor de Cox por pseudo-componente; usa `diffusion_coeff_in_air()` (temperatura en K) y `transfer_coeff_air()`. Con `processes:slick_thickness_model=field` el espesor viene del campo de mancha. Mueve masa de `comp_*` → `comp_*_evaporated` (masa evaporada acumulada por componente) y de `mass_oil` → `mass_evaporated`. |
| `spreading_cicese` | `(self, gravity=9.81, drag_coef=0.003, film_limit=0.000005)` | Esparcimiento gravedad-viscoso (Nihoul, 1984); actualiza `spillet_area`/`spillet_thickness`. `film_limit` es el espesor mínimo de película (m). |
| `emulsification_cicese` | `(self)` | Emulsificación agua-en-crudo tipo Bullwinkle; actualiza `max_water`, `water_fraction`, `viscosity`, `IFT`. Requiere `max_water>0` al sembrar para crudos sin datos de emulsificación en ADIOS (p. ej. MAYA). La tasa efectiva de captación de agua es: `k_emul_eff = k_emul * clip(nu_dry/EMUL_STABILITY_NU_REF, EMUL_STABILITY_MIN, EMUL_STABILITY_MAX) * EMUL_RATE_SCALE`, donde `nu_dry` [cSt] es la viscosidad cinemática del crudo seco (corregida por evaporación) calculada en cada paso. Constantes de clase: `EMUL_STABILITY_NU_REF=1000.0`, `EMUL_STABILITY_MIN=1.0`, `EMUL_STABILITY_MAX=1.0`, `EMUL_RATE_SCALE=1.0` (desde el 2026-10-29 la tasa de Bullwinkle original de NOAA; el ajuste anterior 0.3 / 3.0 / 0.0142 a un objetivo de bolas de alquitrán a 17.2 días no tenía respaldo en datos y se retiró). Véase Manual §8. |
| `dispersion_cicese` | `(self)` | Sumidero de masa de dispersión natural de Roy, que actúa solo sobre elementos superficiales (`z == 0`), condicionado por `_roy_sink_enabled()` (`processes:dispersion` y `processes:dispersion_mode != 'transport'`). Tasa de la función pura `roy_dispersed_fraction()`: `C_Roy · c_disp · V_entrain · dt / (emulsion_density · thickness)`, acotada a 1, con la forma de PyGNOME `sqrt(viscosidad cinemática)`; el espesor es el `oil_film_thickness` del elemento. Mueve masa a `mass_dispersed`. |
| `local_sea_water_density` | `(self)` | Densidad del agua de mar a partir de la temperatura (K → °C) y la salinidad propias de los elementos; se usa en la densidad de la emulsión, en la tasa de arrastre de Li et al. (2017) y en el espectro de gotas (la llamada simple a `sea_water_density()` devuelve el valor por defecto de OpenDrift, 1026.95 kg/m³ a 10 °C, S 35). |
| `surface_wave_mixing` | `(self)` | Sobrescritura: omite el arrastre de Li en `processes:dispersion_mode='budget'`; en otro caso delega a `OpenOil`. |
| `_slick_field_step` | `(self)` | Un paso del espesor de mancha en malla (`processes:slick_thickness_model=field`): reescalar por B/A, advectar con la deriva superficial, esparcir (régimen II de Fay), agregar el volumen recién llegado a la superficie como lentes de Barenblatt; cada elemento superficial lee el espesor de su celda. |
| `handle_subsea_dissolution` | `(self)` | Contabilidad de masa de fase disuelta para partículas `is_dissolved==1`. En el primer encuentro transfiere `cicese_mass_balance['mass_components']` → `mass_dissolved` (una sola vez, poniendo en cero los componentes) y fija `mass_dissolved_subsea` como total de referencia. Véase Manual §7.3 para el ciclo de vida completo (este es el punto de entrada que consume la rama disuelta de `biodegradation_cicese()`). |
| `handle_gas_particles` | `(self)` | Para partículas `is_bubble==1`: registra `mass_gas` y luego desactiva la partícula (estado `released_gas`); la masa de las burbujas de gas se contabiliza pero no se advecta más. |
| `biodegradation_cicese` | `(self)` | **NUEVO.** Decaimiento de primer orden por componente, escalado con Q10 (T_REF=20 °C, `BIODEG_Q10=3`), usando `K_BIODEG_ALI`/`K_BIODEG_AROM`/`K_BIODEG_RESIDUE`, con `BIODEG_SURFACE_FACTOR=1/30` para partículas `z==0` y un factor de escalamiento por tamaño de gota `(BIODEG_DROPLET_D_REF/d)^BIODEG_DROPLET_EXPONENT` (auditoría #2; `d=elements.diameter` acotado a `[BIODEG_DROPLET_DIAMETER_MIN, BIODEG_DROPLET_DIAMETER_MAX]`, `diameter==0`→factor=1), más una rama de fase disuelta (`K_BIODEG_WATER=0.0578 d⁻¹`) para partículas `is_dissolved==1`, que se desactivan cuando `mass_dissolved < BIODEG_WATER_DEPLETION_FRAC` (0.01) de su referencia inicial. Escribe `mass_biodegraded_from_oil/_water`, `mass_biodegraded`, `fraction_biodegraded_from_oil/_water`. Física completa en Manual §7. |
| `photooxidation_cicese` | `(self)` | **NUEVO.** Decaimiento de primer orden por UV de los grupos aromáticos G1-G8 únicamente, usando `K_PHOTO_AROM`, atenuación de Beer-Lambert con la profundidad (`K_D_UV=0.10 m⁻¹`, corte `UV_DEPTH_CUTOFF=50 m`) y `_compute_uv_irradiance()`. De la masa fotooxidada en cada paso, una fracción `OP_DISSOLVED_FRACTION=0.5` se dirige al reservorio disuelto de fotoproductos oxigenados `mass_op_dissolved` en vez de `mass_photooxidized` (auditoría #4); `mass_op_dissolved` decae incondicionalmente (incluso de noche) a `K_OP_REMOVAL=0.116 d⁻¹` (escalada con Q10) hacia `mass_op_degraded`. Escribe `mass_photooxidized`, `mass_op_dissolved`, `mass_op_degraded`, `fraction_photooxidized`, `fraction_op_dissolved`, `fraction_op_degraded`. Física completa en Manual §6. |

### 4.4 Arrastre por oleaje y distribuciones de tamaño de gota

| Método | Firma | Descripción |
|---|---|---|
| `oil_wave_entrainment_rate` | `(self)` | Tasa de arrastre de crudo a la columna de agua por rotura de oleaje (clave `wave_entrainment:entrainment_rate`, actualmente siempre Li et al. 2017). |
| `get_wave_breaking_droplet_diameter_johansen2015` | `(self)` | Espectro de tamaños de gota tras el arrastre por rotura de oleaje, parametrización de Johansen et al. (2015). |
| `get_wave_breaking_droplet_diameter_liz2017` | `(self)` | Espectro de tamaños de gota, parametrización de Li et al. (2017). Se elige con `wave_entrainment:droplet_size_distribution`. |

### 4.5 Dispersantes

| Método | Firma | Descripción |
|---|---|---|
| `dispersant_efficiency_sampled` | `(self, target_DOR, efficiency, samples)` | Muestreo Monte Carlo de la razón dispersante-crudo efectiva dado un DOR objetivo y una distribución de eficiencia. |
| `treated_IFT_and_WC` | `(self, dispersant, applied_DORs)` | Calcula la tensión interfacial y el contenido máximo de agua posteriores al tratamiento para un `dispersant` dado (uno de `OpenCiceseOil.dispersants`: `C9500`, `OSR-52`, `Dasic-NS`, `Average`) y un arreglo de DOR aplicados. |

### 4.6 Auxiliares de propiedades físicas

| Método | Firma | Descripción |
|---|---|---|
| `oil_density` | `(self, indices)` | Densidad del crudo por partícula en las posiciones `indices` del arreglo, a partir de la composición actual. Las filas de `cicese_mass_balance` se buscan por ID de elemento (`self.elements.ID[indices]`; la indexación posicional era un error corregido en la rama `fix/gas-bubble-mask`, manual §10C). |
| `component_molfrac` | `(self, surface)` | Fracciones molares de los 17 pseudo-componentes para los elementos `surface` (máscara booleana); las usa la evaporación. |
| `diffusion_coeff_in_air` | `(self, atm_pressure=101325.)` | Coeficientes de difusión de vapor en aire (Hilal & Karickhoff 2003; Kim & Monroe 2014), por pseudo-componente, función de `atm_pressure`. |
| `transfer_coeff_air` | `(self, surface_indices)` | Coeficientes de transferencia de masa en la interfaz aire-mar para los elementos `surface_indices` (dependen de la rapidez del viento); los usa `evaporation_cicese()`. |
| `_compute_uv_irradiance` | `(self)` | Irradiancia UV solar en superficie (z=0) (W/m²). Prefiere datos reales de lector `surface_downward_uv_radiation` (ERA5), con respaldo en un modelo sintético de coseno diurno. Devuelve 0 de noche. Véase Manual §6.2. |

### 4.7 Salida, balances y análisis

| Método | Firma | Descripción |
|---|---|---|
| `get_oil_budget` | `(self)` | Devuelve un diccionario con `oil_density`, `mass_dispersed`, `mass_submerged`, `mass_surface`, `mass_stranded`, `mass_dissolved_subsea`, `mass_gas`, `mass_evaporated`, `mass_biodegraded`, `mass_photooxidized`, `mass_total`, cada uno una serie temporal 1-D (sumada sobre las trayectorias). `mass_total` debe igualar la masa liberada (verificación de conservación de masa). Devuelve `None` en simulaciones hacia atrás. Usa el auxiliar `_prop()` para compatibilidad de xarray con OD 1.14.9. |
| `plot_oil_budget` | `(self, filename=None, ax=None, show_density_viscosity=True, show_wind_and_current=True)` | Gráfica de áreas apiladas de las categorías de `get_oil_budget()` en el tiempo, opcionalmente con subpaneles de densidad/viscosidad y viento/corriente. Guarda en `filename` o muestra de forma interactiva. |
| `plot_density` | `(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None, cmapname=...)` | Grafica un campo de densidad en malla (p. ej. de `get_density_array_framed()`) en un mapa. |
| `plot_stats_property` | `(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None, cmapname=...)` | Grafica una propiedad estadística en malla (p. ej. de `get_property_binned()`) en un mapa. |
| `create_new_map` | `(self, corners=None, buffer=0.1, lscale='auto', fast=False)` | Crea una `Figure`/`GeoAxes` de Cartopy para gráficas de trayectorias/densidad, con rasgos de tierra a la resolución `lscale`. |
| `get_bins` | `(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, depths=None)` | Construye una malla regular de lon/lat (y opcionalmente de profundidad) de tamaño `bin_size_deg` para agregación espacial. |
| `get_property_binned` | `(self, bins, p_groups, depths, de_sizes, property='particles_no')` | Agrega una `property` de elemento con nombre (por defecto: conteo de partículas) en la malla `bins`, opcionalmente separada por grupo de pseudo-componente (`p_groups`), bin de profundidad y bin de tamaño de gota. |
| `get_unique_visits` | `(self, bins, p_groups, depths, de_sizes)` | Cuenta visitas únicas de partículas por bin espacial (evita contar doble a una partícula que revisita una celda). |
| `get_density_array_framed` | `(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, weight=None)` | Envoltura de conveniencia que combina `get_bins()` + `get_property_binned()` para producir un arreglo 2-D de densidad de masa/partículas en un marco fijo de lon/lat, opcionalmente ponderado por `weight` (p. ej. `mass_oil`). |

---

## 5. `fluid_properties.FluidProps`

`class FluidProps(object)` en `fluid_properties.py` — el marco de
pseudo-componentes de 17 ("cortes de destilación") que sustenta la
intemperización por componentes de CICOIL (evaporación, biodegradación,
fotooxidación).

| Miembro | Firma | Descripción |
|---|---|---|
| `_riazi` | `(T, A, B)` (función de módulo) | Auxiliar de la correlación de Riazi: ajuste genérico tipo `Y = exp(A + B*T)` que usan varios de los estimadores de propiedades siguientes. |
| `__init__` | `(self, oiltype, max_cuts)` | Construye el marco de componentes para `oiltype` (un objeto `Oil`/`CiceseOil` de ADIOS) con hasta `max_cuts` cortes de destilación (CICOIL usa 17 = 8 alifáticos + 8 aromáticos + 1 residuo). |
| `_create_TBP` | `(self, max_cuts)` | Construye la curva de destilación de punto de ebullición verdadero (TBP) y los puntos de ebullición por corte, usados para derivar pesos moleculares, densidades y presiones de vapor de cada pseudo-componente. |
| `get_Riazi_model` | `(self)` | Devuelve el conjunto de propiedades basado en la correlación de Riazi (peso molecular, densidad, propiedades críticas) de cada pseudo-componente. |
| `get_live_composition` | `(self, GOR=0)` | Devuelve la composición "viva" (in situ, ajustada por presión) considerando la razón gas-crudo `GOR`; relevante para liberaciones submarinas/TAMOC donde el gas disuelto afecta las propiedades globales. |
| `setup_treatment` | `(self, treatment, dispersant, DOR, disp_stock, delay, time_frac_surf, time_frac_ssdi, effic_surf, effic_ssdi, wind_limit, equipment_rate)` | Configura un escenario de tratamiento con dispersante químico (usado por `seed_with_dispersant()`/`OpenCiceseOil`); guarda los parámetros del tratamiento para posteriores ajustes de IFT/max-water. |
| `verify_cut_fractional_masses` | `(cls, fmass_i, T_i, f_sat_i, f_arom_i, prev_f_sat_i)` (método de clase) | Verifica que las masas fraccionales por corte (`fmass_i`), los puntos de ebullición (`T_i`) y las particiones saturados/aromáticos (`f_sat_i`, `f_arom_i`) sean físicamente consistentes (monótonas, suman 1, etc.) dada la fracción de saturados del corte anterior. |

---

## 6. `cicoil_estimations.py` — herramientas de estimación para vacíos de ADIOS

Reemplazo independiente de `adios.util.estimations` (ausente en `adios_db`). 32
funciones, todas puras/sin estado, en su mayoría correlaciones empíricas de la
documentación de ADIOS2, Riazi, Fingas o las recomendaciones de Bill Lehr.
Agrupadas por propósito:

### 6.1 Densidad / gravedad API / gravedad específica
| Función | Firma |
|---|---|
| `density_from_api` | `(api)` |
| `api_from_density` | `(density)` |
| `density_at_temp` | `(ref_density, ref_temp_k, temp_k, k_rho_t=0.0008)` |
| `vol_expansion_coeff` | `(rho_0, t_0, rho_1, t_1)` |
| `specific_gravity` | `(density)` |

### 6.2 Viscosidad
| Función | Firma |
|---|---|
| `dvis_to_kvis` | `(dvis, density)` — viscosidad dinámica → cinemática |
| `kvis_at_temp` | `(ref_kvis, ref_temp_k, temp_k, k_v2=2416.0)` |

### 6.3 Fracciones SARA (saturados/aromáticos/resinas/asfaltenos)
| Función | Firma |
|---|---|
| `resin_fraction` | `(density, viscosity, f_other=0.0)` |
| `asphaltene_fraction` | `(density, viscosity, f_other=0.0)` |
| `saturates_fraction` | `(density, viscosity, f_other=0.0)` |
| `aromatics_fraction` | `(f_res, f_asph, f_sat)` |
| `_A_coeff` | `(density)` — coeficiente empírico de Fingas |
| `_B_coeff` | `(density, viscosity)` — coeficiente empírico de Fingas |

### 6.4 Cortes de destilación
| Función | Firma |
|---|---|
| `cut_temps_from_api` | `(api, N=5)` — puntos de ebullición de N cortes a partir de la gravedad API |
| `fmasses_from_cuts` | `(f_evap_i)` — masas fraccionales a partir de los cortes de la curva de evaporación |
| `fmasses_flat_dist` | `(f_res, f_asph, N=5)` — respaldo de distribución plana para N cortes |

### 6.5 Peso molecular y densidad por corte (correlaciones de Riazi)
| Función | Firma |
|---|---|
| `saturate_mol_wt` | `(boiling_point)` |
| `aromatic_mol_wt` | `(boiling_point)` |
| `resin_mol_wt` | `(_boiling_point)` |
| `asphaltene_mol_wt` | `(_boiling_point)` |
| `trial_densities` | `(boiling_points, watson_factor)` |
| `saturate_densities` | `(boiling_points)` |
| `aromatic_densities` | `(boiling_points)` |
| `resin_densities` | `(_boiling_points)` |
| `asphaltene_densities` | `(_boiling_points)` |

### 6.6 Caracterización de hidrocarburos
| Función | Firma |
|---|---|
| `_hydrocarbon_characterization_param` | `(specific_gravity, temp_k)` |
| `refractive_index` | `(hc_char_param)` |
| `_hydrocarbon_grouping_param` | `(mol_wt, specific_gravity, temp_k)` |
| `saturate_mass_fraction` | `(fmass_i, mol_wt, specific_gravity, temp_k)` |

### 6.7 Tensión interfacial, punto de fluidez, punto de inflamación, Bullwinkle
| Función | Firma |
|---|---|
| `oil_water_surface_tension_from_api` | `(api)` |
| `pour_point_from_kvis` | `(ref_kvis, ref_temp_k)` |
| `pour_point_from_sg_mw_kvis` | `(specific_gravity, mol_wt, kvis)` |
| `flash_point_from_bp` | `(temp_k)` |
| `flash_point_from_api` | `(api)` |
| `bullwinkle_fraction_from_asph` | `(f_asph)` — fracción de inicio de emulsificación a partir del contenido de asfaltenos |
| `bullwinkle_fraction_from_api` | `(api)` — fracción de inicio de emulsificación a partir de la gravedad API |

---

## 7. Acoplamiento con TAMOC

### 7.1 `tamoc_plume.Plume`
`class Plume(object)` en `tamoc_plume.py` — representa la salida de una corrida
del modelo de pluma curva (BPM) de TAMOC, empaquetada para la siembra submarina
en OpenDrift.

| Método | Firma | Descripción |
|---|---|---|
| `__init__` | `(self, chem_data_file=None, plume_file=None)` | Opcionalmente carga de inmediato los archivos de datos químicos y de salida de la pluma. |
| `load_chem_data` | `(self, fname=None)` | Carga las tablas de propiedades químicas de TAMOC (véase `tamoc_chemical_properties.py`). |
| `load_plume` | `(self, fname=None)` | Carga la salida de trayectoria de la pluma del BPM de TAMOC (posición de la línea central, radio, velocidad, temperatura, fracciones de fase frente a la altura sobre la liberación). |
| `get_particle_properties` | `(self, particles_file, tracked=False)` | Lee la salida por partícula de TAMOC (`particles_file`) y devuelve arreglos de posiciones, tamaños, fase (gota de crudo / burbuja de gas / disuelto) y masa; los consume directamente `OpenCiceseOil.seed_plume_elements()`. `tracked=True` conserva las trayectorias de partículas individuales en lugar de agregarlas. |

### 7.2 `run_tamoc.run_tamoc`
```python
run_tamoc(live_comp, composition, chemdata, chemunits, output_name,
          profile, total_flow, z0, D, Tj, phi_0, theta_0, bins)
```
Envoltura que corre una simulación BPM de TAMOC dados:
- `live_comp` / `composition` / `chemdata` / `chemunits` — composición del fluido y unidades (de `FluidProps.get_live_composition()` y `tamoc_chemical_properties.py`)
- `profile` — perfil oceánico ambiental (T, S, corrientes frente a la profundidad)
- `total_flow`, `z0`, `D`, `Tj`, `phi_0`, `theta_0` — tasa de liberación, profundidad, diámetro del orificio, temperatura del chorro y ángulos de descarga
- `bins` — número de bins de tamaño de partícula de la fase dispersa
- `output_name` — nombre base de los archivos de salida de TAMOC, que luego leen `Plume.load_plume()` / `Plume.get_particle_properties()`

### 7.2b `vertical_profile` (rama `feat/vertical-profile`)
```python
sample_profile(readers, lon, lat, time, z_max=None, add_surface=True)
create_profile(nc_file, readers, lon, lat, time, source='OpenDrift readers',
               summary='Ambient profile from OpenDrift readers', sea_name=..., z_max=None, add_surface=True)
```
`sample_profile` lee T, S, u, v en cada nivel de profundidad de los lectores (columna cruda de la celda más cercana con `reader.get_variables`; se conservan los NaN bajo el fondo) y devuelve las columnas. `create_profile` escribe el netCDF de perfil ambiental de TAMOC (T, S, presión, ua, va contra profundidad) y devuelve su ruta; páselo como `profile` a `run_tamoc`. Los lectores necesitan niveles verticales (archivos CROCO acoplados; los GLORYS/ANFC de Abkatun 2026 tienen uno). Use `Plume(tracked=False)` con TAMOC 4.1.2.

### 7.3 `tamoc_chemical_properties.py`
Tablas estáticas de propiedades (pesos moleculares, propiedades críticas,
constantes de Henry, etc.) de los componentes químicos que usan los cálculos de
equilibrio y disolución de TAMOC; alimenta a `run_tamoc()`.

---

## 8. Lectores NEMO/CROCO

Ambos lectores activos heredan de `opendrift.readers.basereader.BaseReader` y
exponen los campos oceánicos mediante un diccionario `NEMO_variable_mapping`
(variable nativa → `standard_name` de CF), que la interfaz de lectores de
OpenDrift consume automáticamente.

### 8.1 `reader_NEMO_native_v3.Reader` (activo — nombres nativos de CROCO)
```python
Reader(filename=None, filenameU=None, filenameV=None, name=None,
       gridfile=None, custom_var_mapping=None, _FillValue=None)
```
Mapeo de variables (nombres diagnósticos nativos de CROCO/NEMO):

| Nativo | `standard_name` de CF |
|---|---|
| `sossheig` | `sea_surface_height` |
| `vozocrtx` | `x_sea_water_velocity` |
| `vomecrty` | `y_sea_water_velocity` |
| `vovecrtz` | `upward_sea_water_velocity` |
| `votemper` | `sea_water_temperature` |
| `vosaline` | `sea_water_salinity` |
| `utau` | `surface_downward_x_stress` |
| `vtau` | `surface_downward_y_stress` |

`custom_var_mapping` (dict) sobrescribe/extiende este mapeo.

### 8.2 `reader_nemo_optimized.Reader` (activo — nombres genéricos NEMO/CMEMS)
```python
Reader(filename=None, filenameU=None, filenameV=None, filenameW=None,
       filename_mask=None, name=None, gridfile=None,
       custom_var_mapping=None, _FillValue=None)
```
Mapeo de variables (nombres cortos estilo CMEMS):

| Nativo | `standard_name` de CF |
|---|---|
| `ssh` | `sea_surface_height` |
| `uoce` | `x_sea_water_velocity` |
| `voce` | `y_sea_water_velocity` |
| `wo` | `upward_sea_water_velocity` |
| `toce` | `sea_water_temperature` |
| `soce` | `sea_water_salinity` |
| `tauuo` | `surface_downward_x_stress` |
| `tauvo` | `surface_downward_y_stress` |
| `tmask` | `land_binary_mask` |

Admite archivos separados de U/V/W/máscara (`filenameU/V/W`, `filename_mask`)
además de un `filename` combinado.

### 8.3 `reader_nemo_combined.py` / `reader_nemo_modified.py` — lectores NEMO alternativos
Ambos archivos se limpiaron (2026-06-14) quitando un borrador anterior
totalmente comentado que precedía a la `class Reader` funcional de cada
archivo. Las clases vivas restantes usan el mismo mapeo
`ssh`/`uoce`/`voce`/`wo`/`toce`/`soce`/`tauuo`/`tauvo`/`tmask` que
`reader_nemo_optimized.py`, con soporte de archivos separados U/V/W/máscara
(`filenameU/V/W`, `filename_mask`) y un `gridfile` opcional para
`nav_lon`/`nav_lat`:

- `reader_nemo_combined.py` — basado en xarray (`xr.open_dataset`, carga los
  conjuntos de datos en memoria, interpola con `scipy.ndimage.map_coordinates`).
  Se importa y se ejercita en `test_cicoil_e2e.py` (T5.3).
- `reader_nemo_modified.py` — basado en netCDF4 (`Dataset`/`MFDataset`, admite
  patrones glob mediante `MFDataset`), con un auxiliar `lonlat_to_ij()` para
  buscar el punto de malla más cercano. Actualmente no lo usan la batería de
  pruebas ni los scripts de producción; se conserva como implementación
  alternativa basada en netCDF4.

`reader_nemo_optimized.py` y `reader_NEMO_native_v3.py` siguen siendo los
lectores principales para corridas de producción.

### 8.4 `reader_croco_native.Reader` (activo — salida cruda de CROCO en coordenada sigma)
```python
Reader(filename, name=None, gridfile=None, standard_name_mapping={}, destagger=True, **kwargs)
```
Subclase de `opendrift.readers.reader_ROMS_native.Reader` (OpenDrift 1.14.9),
por lo que se aplican todos los argumentos y mapeos de variables de
`reader_ROMS_native`. Lee `croco_avg.nc` / `croco_his.nc` directamente (niveles
sigma, malla C escalonada, `Cs_rho`, `hc`, `Vtransform`, `angle`,
`lon_rho`/`lat_rho`) y convierte sigma a profundidad al vuelo.

| Miembro | Descripción |
|---|---|
| `destagger` (argumento, por defecto `True`) | Promedia `u` (xi_u) y `v` (eta_v) a los puntos rho de forma perezosa (dask) tras abrir, copiando los puntos extremos; `False` reproduce el lector original, que los usa como si estuvieran en los nodos rho (desplazamiento de media celda, ~2.3 km en una malla de 1/24 de grado). |
| `mask_u`, `mask_v` (propiedades) | Derivadas de `mask_rho` (húmedo solo si ambos puntos rho vecinos son húmedos); la salida de CROCO no las trae y el lector original falla con `KeyError: 'mask_u'`. |
| `destaggered` | `True` una vez aplicado el desescalonamiento. |

Diferencias con `reader_NEMO_native_v3` (§8.1): ese lector espera archivos de
niveles z estilo NEMO con nombres de variable de NEMO y no maneja coordenadas
sigma. Solo variables oceánicas; el viento y el oleaje requieren otros
lectores. El primer registro de `croco_avg.nc` está a la mitad de la primera
ventana de promedio, por lo que una corrida no puede empezar antes. Validado
para transporte superficial (véase Manual §3 y §9); los campos subsuperficiales
aún no están validados. Se registra una advertencia cuando OpenDrift no es
1.14.9. Prueba: `test_croco_native_reader.py`.

---

## 9. `export/io_stat_nc.py` — salida estadística en malla

Función a nivel de módulo (se parcha en una instancia del modelo OpenDrift como
gancho de exportación, convención `io_*` de NOAA/OpenDrift):

| Función | Firma | Descripción |
|---|---|---|
| `init` | `(self, model_template, filename, prop_metadata)` | Abre un nuevo archivo NetCDF (`filename`) para salida estadística en malla, escribiendo los atributos globales (`model_url`, `opendrift_class`, `opendrift_module`, `readers`, `time_step_calculation`) y las definiciones de dimensiones/variables derivadas de `prop_metadata` (las propiedades en malla que se escriben en cada paso de tiempo, p. ej. de `get_density_array_framed()` / `get_property_binned()`). |

Se usa para producir campos resumen en malla (con bins de lon/lat/profundidad)
en lugar de (o además de) la salida NetCDF por trayectoria; útil en ensambles
grandes donde la salida por partícula sería demasiado grande.

---

## 10. `plotting.py`

8 auxiliares de graficación de calidad de publicación, todos funciones
independientes (operan sobre diccionarios de `get_oil_budget()`, archivos de
salida NetCDF o una instancia de `OpenCiceseOil`; sin estado compartido).

| Función | Firma | Descripción |
|---|---|---|
| `plot_mass_budget` | `(budget, title=None, timestep_hours=1, unit='tonnes', show_photoox=True, filename=None, figsize=..., dpi=...)` | Gráfica de áreas apiladas de un diccionario de balance de crudo (de `get_oil_budget()` o de un `.npz`/`.csv` guardado) en el tiempo, en `unit` (toneladas/kg/%). `show_photoox` activa la categoría fotooxidada. |
| `plot_trajectory_map` | `(nc_file, spill_lon=None, spill_lat=None, title=None, time_index=-1, color_by=None, filename=None, figsize=..., dpi=...)` | Mapa de Cartopy de las trayectorias de partículas de un `nc_file` de salida de OpenDrift, opcionalmente coloreado por `color_by` (p. ej. `'status'`, `'z'`, una variable `comp_*`), marcando el origen del derrame. |
| `plot_photoox_timeseries` | `(budget_with, budget_without, timestep_hours=1, title=None, filename=None, figsize=..., dpi=...)` | Compara las series de tiempo de `fraction_photooxidized` entre una corrida con fotooxidación por UV (`budget_with`) y una corrida de control sin ella (`budget_without`). |
| `plot_aromatic_budget` | `(o, title=None, filename=None, figsize=..., dpi=...)` | Gráfica de barras apiladas de las fracciones de masa evaporada + fotooxidada por grupo aromático (G1–G8), leyendo `comp_arom_G*` / `comp_arom_G*_evaporated` / totales de fotooxidación directamente de una instancia `o` de `OpenCiceseOil`. |
| `plot_f_surface_sensitivity` | `(Hs_range, Tp_range, w_rise=..., n=..., title=None, filename=None, figsize=..., dpi=...)` | Mapa de calor de la corrección `f_surface` de sumergimiento por oleaje (antigua, ya eliminada del modelo) en función de la altura significativa `Hs_range` y el periodo pico `Tp_range`; se conserva para comparación diagnóstica/histórica (véase Manual §6.2 sobre por qué se eliminó `f_surface` de `photooxidation_cicese()`). |
| `plot_surface_density` | `(nc_file, time_index, spill_lon=None, spill_lat=None, bin_deg=..., weight_var='mass_oil', title=None, filename=None, figsize=..., dpi=...)` | Mapa de densidad en malla de la masa superficial de crudo (u otra `weight_var`) en el paso de tiempo `time_index`, con bins de resolución `bin_deg`. |
| `compare_weathering_modes` | `(budgets, labels, timestep_hours=1, title=None, filename=None, figsize=..., dpi=...)` | Comparación de varios paneles lado a lado de las fracciones del destino del crudo entre varios `budgets` (p. ej. con/sin biodegradación, con/sin fotooxidación), cada uno etiquetado según `labels`. |
| `plot_photoox_rate_constants` | `(filename=None, figsize=..., dpi=...)` | Gráfica de barras de las constantes de tasa / vidas medias actuales de `K_PHOTO_AROM` G1–G8 (lee directamente de `ciceseoil.OpenCiceseOil.K_PHOTO_AROM`); útil para documentar la calibración actual (Manual §6.3). |

---

## 11. `download_data.py` — adquisición y preprocesamiento de forzamientos

| Función | Firma | Descripción |
|---|---|---|
| `download_glorys12` | `(lon_min, lon_max, lat_min, lat_max, date_start, date_end, outfile, depth_min=..., depth_max=...)` | Descarga corrientes/T/S oceánicas de la reanálisis GLORYS12 (Copernicus Marine) para la caja, el rango de fechas y el rango de profundidad dados; recomendado para corridas de hindcast reproducibles (véase la advertencia del Manual §9 sobre ANFC). |
| `download_mercator_anfc` | `(lon_min, lon_max, lat_min, lat_max, date_start, date_end, outfile, depth_min=..., depth_max=..., tmpdir=...)` | Descarga campos oceánicos de análisis/pronóstico casi en tiempo real de Mercator ANFC; úsese solo cuando no se requiere reproducibilidad entre fechas de descarga. |
| `download_era5_winds` | `(lon_min, lon_max, lat_min, lat_max, year, months, outfile, time_step=...)` | Descarga las componentes del viento a 10 m de ERA5 (y campos superficiales relacionados) mediante la API de CDS. |
| `download_era5_waves` | `(lon_min, lon_max, lat_min, lat_max, year, months, outfile, time_step=...)` | Descarga campos de oleaje de ERA5 (Hs, Tp, deriva de Stokes) para cálculos de sumergimiento por oleaje y dispersión. |
| `_preprocess_era5_cf` | `(ds)` | Interna: agrega atributos `standard_name` de CF y normaliza la dimensión/unidades de tiempo en un `xarray.Dataset` crudo de ERA5 para que `reader_netCDF_CF_generic` de OpenDrift reconozca las variables (necesario: el NetCDF crudo de ERA5 trae `standard_name='unknown'`). |
| `preprocess_era5_zip` | `(zip_nc_path, outfile)` | Desempaca un archivo NetCDF de ERA5 disfrazado de ZIP (formato de entrega antiguo de la API de CDS) y aplica `_preprocess_era5_cf`. |
| `check_domain_coverage` | `(forcing_files, spill_lon, spill_lat, sim_days, buffer_deg=...)` | Verifica que los `forcing_files` dados cubran, en espacio y tiempo, una simulación que empieza en `(spill_lon, spill_lat)` durante `sim_days`, con un margen `buffer_deg`; úsese antes de lanzar una simulación larga para detectar temprano fallas por falta de cobertura. |
| `list_croco_files` | `(croco_dir, pattern=...)` | Lista los archivos NetCDF de salida de CROCO en `croco_dir` que coinciden con `pattern`, para usarlos con `reader_NEMO_native_v3.Reader`. |

> Para la radiación UV de ERA5 en particular (forzamiento de fotooxidación),
> véase `preprocess_era5_uv.py` (no es parte de `download_data.py`): produce los
> archivos `surface_downward_uv_radiation` conformes a CF en
> `data_local/era5_uv_cf/` a los que se hace referencia en el Manual §6.2.

---

**Versión:** CICOILv2 (commit `b032140`; etiqueta de versión `v2.1.0` pendiente)
**Última actualización:** 2026-11-17
