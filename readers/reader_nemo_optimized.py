#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reader NEMO completo con proyección - Versión final corregida
"""
from opendrift.readers.basereader import BaseReader
from netCDF4 import Dataset, num2date
from pyproj import CRS
import numpy as np
import logging
import os
import time  # Solo para time.time()

class Reader(BaseReader):
    def __init__(self, filename=None, filenameU=None, filenameV=None,
                 filenameW=None, filename_mask=None, name=None,
                 gridfile=None, custom_var_mapping=None, _FillValue=None):

        if filename is None:
            raise ValueError('Debes proporcionar al menos el archivo de grid T')

        # Configuración de variables
        self.NEMO_variable_mapping = {
            'ssh': 'sea_surface_height',
            'uoce': 'x_sea_water_velocity',
            'voce': 'y_sea_water_velocity',
            'wo': 'upward_sea_water_velocity',
            'toce': 'sea_water_temperature',
            'soce': 'sea_water_salinity',
            'tauuo': 'surface_downward_x_stress',
            'tauvo': 'surface_downward_y_stress',
            'tmask': 'land_binary_mask'
        }
        if custom_var_mapping:
            self.NEMO_variable_mapping.update(custom_var_mapping)

        self.name = name or "NEMO_completo"
        self._FillValue = _FillValue

        # Solo guardamos las rutas
        self.files = {
            'T': filename,
            'U': filenameU,
            'V': filenameV,
            'W': filenameW,
            'mask': filename_mask
        }
        
        # Abrir solo para metadatos
        try:
            start_meta = time.time()  # Usar time.time() para timestamp
            with Dataset(filename, 'r') as ds:
                # Dimensiones básicas
                self.lat = ds.variables['nav_lat'][:]
                self.lon = ds.variables['nav_lon'][:]
                
                # Tiempos
                ocean_time = ds.variables['time_counter']
                time_units = ocean_time.units
                self.times = num2date(ocean_time[:], time_units)
                self.start_time = self.times[0]
                self.end_time = self.times[-1]
                
                # Profundidades
                if 'deptht' in ds.variables:
                    varz = ds.variables['deptht']
                    positive_attr = getattr(varz, 'positive', 'down')
                    self.z = -varz[:] if positive_attr != 'up' else varz[:]
                else:
                    self.z = None
            
            logging.info(f"Metadatos cargados en {time.time() - start_meta:.2f} segundos")
        except Exception as e:
            raise ValueError(f"Error leyendo metadatos: {e}")

        # Configuración de la rejilla
        self.xmin = np.min(self.lon)
        self.xmax = np.max(self.lon)
        self.ymin = np.min(self.lat)
        self.ymax = np.max(self.lat)
        self.delta_x = np.mean(np.abs(np.diff(self.lon[0])))
        self.delta_y = np.mean(np.abs(np.diff(self.lat[:,0])))
        
        # Configuración CRS crítica - WGS84 (EPSG:4326)
        self.proj4 = '+proj=latlong +datum=WGS84'
        self.crs = CRS.from_proj4(self.proj4)
        
        self.variables = list(self.NEMO_variable_mapping.values())
        
        super().__init__()
        logging.info(f"Lector {self.name} inicializado. Proyección: {self.proj4}")
        logging.info(f"Dominio: {self.xmin:.2f}-{self.xmax:.2f}E, {self.ymin:.2f}-{self.ymax:.2f}N")

    def _get_dataset(self, grid_type):
        """Abre el dataset cuando sea necesario"""
        file_path = self.files.get(grid_type)
        if not file_path or not os.path.exists(file_path):
            return None
        
        try:
            return Dataset(file_path, 'r')
        except Exception as e:
            logging.error(f"Error abriendo {file_path}: {e}")
            return None

    def _get_grid_type_for_variable(self, var_name):
        """Determina el tipo de grid para una variable"""
        if var_name in ['uoce', 'tauuo']:
            return 'U'
        elif var_name in ['voce', 'tauvo']:
            return 'V'
        elif var_name in ['wo']:
            return 'W'
        elif var_name in ['tmask']:
            return 'mask'
        else:
            return 'T'  # Grid T por defecto

    def _read_variable(self, ds, varname, time_index, depth_index=0):
        """Lee una variable específica con manejo de dimensiones"""
        if varname not in ds.variables:
            logging.warning(f"Variable {varname} no encontrada en dataset")
            return None
        
        var = ds.variables[varname]
        dims = var.dimensions
        
        # Determinar índices para cada dimensión
        slices = []
        for dim in dims:
            if dim == 'time_counter':
                slices.append(time_index)
            elif dim in ['deptht', 'depthu', 'depthv', 'depthw']:
                slices.append(depth_index)
            else:
                slices.append(slice(None))
        
        try:
            data = var[tuple(slices)]
            # Manejar posibles escalas
            if hasattr(var, 'scale_factor'):
                data = data * var.scale_factor
            if hasattr(var, 'add_offset'):
                data = data + var.add_offset
            return data
        except IndexError as e:
            logging.error(f"Error de índice en {varname}: {e}")
            return None
        except Exception as e:
            logging.error(f"Error leyendo {varname}: {str(e)}")
            return None

    def get_variables(self, requested_variables, time=None, x=None, y=None, z=None, block=False):
        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)
        
        nearestTime, tidx, _, _, _, _ = self.nearest_time(time)
        variables = {}
        
        # Si z no se proporciona, usar 0 (superficie)
        if z is None:
            z = np.zeros(len(x))
        
        # Pre-cálculo de índices espaciales
        logging.info(f"Calculando índices para {len(x)} puntos...")
        start_idx = time.time()  # Usamos time.time() para obtener timestamp en segundos

        # Convertir coordenadas a índices
        lon_flat = self.lon.ravel()
        lat_flat = self.lat.ravel()
        indices = np.zeros(len(x), dtype=int)
        
        for i, (lon, lat) in enumerate(zip(x, y)):
            # Búsqueda eficiente del punto más cercano
            dist = (lon_flat - lon)**2 + (lat_flat - lat)**2
            idx = np.argmin(dist)
            indices[i] = idx
        
        # Calculamos la duración usando time.time()
        duration_sec = time.time() - start_idx
        logging.info(f"Índices calculados en {duration_sec:.2f} segundos")
        
        # Procesar cada variable
        for var in requested_variables:
            # Encontrar el nombre de la variable NEMO
            varname_nemo = next((k for k, v in self.NEMO_variable_mapping.items() if v == var), None)
            if varname_nemo is None:
                logging.warning(f"Variable {var} no encontrada en el mapeo. Saltando...")
                variables[var] = np.ma.masked_array([np.zeros(len(x))], mask=True)
                continue
            
            # Determinar tipo de grid
            grid_type = self._get_grid_type_for_variable(varname_nemo)
            
            # Obtener dataset
            ds = self._get_dataset(grid_type)
            if ds is None:
                logging.warning(f"Dataset no disponible para {var} ({grid_type})")
                variables[var] = np.ma.masked_array([np.zeros(len(x))], mask=True)
                continue
            
            # Leer datos
            try:
                start_read = time.time()
                
                # Manejo seguro del índice de tiempo
                if 'time_counter' in ds.variables:
                    time_index = min(tidx, len(ds.variables['time_counter'][:]) - 1)
                else:
                    time_index = 0
                
                # Leer datos
                data = self._read_variable(ds, varname_nemo, time_index)
                
                if data is None:
                    logging.warning(f"Variable {varname_nemo} no se pudo leer")
                    variables[var] = np.ma.masked_array([np.zeros(len(x))], mask=True)
                    ds.close()
                    continue
                
                # Aplanar y obtener valores en los índices
                if data.ndim > 1:
                    data_flat = data.ravel()
                else:
                    data_flat = data
                
                # Manejar posibles problemas de tamaño
                if len(data_flat) <= np.max(indices):
                    logging.error(f"Índice fuera de rango: {np.max(indices)} > {len(data_flat)-1}")
                    values = np.full(len(x), np.nan)
                else:
                    values = data_flat[indices]
                
                # Manejar valores faltantes
                fill_value = getattr(ds.variables[varname_nemo], '_FillValue', np.nan)
                if fill_value is not None:
                    values = np.where(values == fill_value, np.nan, values)
                
                # Crear array 2D para OpenDrift [time, particles]
                variables[var] = np.ma.masked_invalid(np.atleast_2d(values))
                
                # Diagnóstico
                non_zero = np.count_nonzero(values)
                read_duration = time.time() - start_read
                logging.info(f"{var}: {non_zero}/{len(values)} no-cero, "
                             f"rango [{np.nanmin(values):.4f}, {np.nanmax(values):.4f}] "
                             f"tiempo lectura: {read_duration:.2f}s")
                
            except Exception as e:
                logging.error(f"Error leyendo {varname_nemo}: {str(e)}")
                variables[var] = np.ma.masked_array([np.zeros(len(x))], mask=True)
            finally:
                if ds:
                    ds.close()

        # Asegurar que todas las variables solicitadas existan
        for var in requested_variables:
            if var not in variables:
                logging.warning(f"Creando array vacío para {var}")
                variables[var] = np.ma.masked_array([np.zeros(len(x))], mask=True)
        
        # Variables requeridas por OpenDrift
        variables['x'] = x
        variables['y'] = y
        variables['z'] = z
        variables['time'] = nearestTime

        return variables