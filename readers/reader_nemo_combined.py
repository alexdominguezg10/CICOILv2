#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reader NEMO combinado personalizado para OpenDrift - Versión corregida
"""
from opendrift.readers.basereader import BaseReader
from netCDF4 import Dataset, num2date
from pyproj import CRS
import numpy as np
import logging
import xarray as xr
import scipy.ndimage

class Reader(BaseReader):
    def __init__(self, filename=None, filenameU=None, filenameV=None,
                 filenameW=None, filename_mask=None, name=None,
                 gridfile=None, custom_var_mapping=None, _FillValue=None):

        if filename is None:
            raise ValueError('Necesitas proporcionar filename (grid T)')

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

        if custom_var_mapping is not None:
            for var_name in list(self.NEMO_variable_mapping.keys()):
                for cvar_name in custom_var_mapping.keys():
                    if custom_var_mapping[cvar_name] == self.NEMO_variable_mapping[var_name]:
                        self.NEMO_variable_mapping.pop(var_name)
                        break
            self.NEMO_variable_mapping.update(custom_var_mapping)

        self.name = name or "NEMO_combined_reader"
        self._FillValue = _FillValue

        # Abrir datasets usando xarray para mejor manejo
        try:
            self.ds = xr.open_dataset(filename, decode_times=False)
            self.dsU = xr.open_dataset(filenameU, decode_times=False) if filenameU else None
            self.dsV = xr.open_dataset(filenameV, decode_times=False) if filenameV else None
            self.dsW = xr.open_dataset(filenameW, decode_times=False) if filenameW else None
            self.ds_mask = xr.open_dataset(filename_mask, decode_times=False) if filename_mask else None
            
            # Cargar datos a memoria
            self.ds.load()
            if self.dsU: self.dsU.load()
            if self.dsV: self.dsV.load()
            if self.dsW: self.dsW.load()
            if self.ds_mask: self.ds_mask.load()

        except Exception as e:
            raise ValueError(f"Error abriendo archivos: {e}")

        # Manejar coordenadas
        if gridfile:
            gf = xr.open_dataset(gridfile)
            self.lat = gf['nav_lat'].values
            self.lon = gf['nav_lon'].values
            gf.close()
        else:
            self.lat = self.ds['nav_lat'].values
            self.lon = self.ds['nav_lon'].values

        # Obtener tiempos
        ocean_time = self.ds['time_counter']
        time_units = ocean_time.attrs.get('units', 'seconds since 1950-01-01')
        self.times = num2date(ocean_time.values, time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.time_step = self.times[1] - self.times[0] if len(self.times) > 1 else None

        # Manejar profundidades
        for dvn in ['deptht', 'depthu', 'depthv', 'depthw']:
            if dvn in self.ds:
                varz = self.ds[dvn]
                positive_attr = varz.attrs.get('positive', 'down')
                self.z = -varz.values if positive_attr != 'up' else varz.values
                break
        else:
            self.z = None

        # Cargar máscaras si existen
        self.tmask = self.ds_mask['tmask'].values if self.ds_mask is not None and 'tmask' in self.ds_mask else None
        self.umask = self.ds_mask['umask'].values if self.ds_mask is not None and 'umask' in self.ds_mask else None
        self.vmask = self.ds_mask['vmask'].values if self.ds_mask is not None and 'vmask' in self.ds_mask else None

        # Configuración de la rejilla
        self.xmin = np.nanmin(self.lon)
        self.xmax = np.nanmax(self.lon)
        self.ymin = np.nanmin(self.lat)
        self.ymax = np.nanmax(self.lat)
        self.delta_x = np.nanmean(np.abs(np.diff(self.lon, axis=1)))
        self.delta_y = np.nanmean(np.abs(np.diff(self.lat, axis=0)))
        
        # Crear rejilla regular para interpolación
        self.xi = np.linspace(self.xmin, self.xmax, self.lon.shape[1])
        self.yi = np.linspace(self.ymin, self.ymax, self.lon.shape[0])
        
        super().__init__()
        self.proj4 = CRS.from_epsg(4326).to_proj4()
        self.variables = list(self.NEMO_variable_mapping.values())

    @property
    def variables_available(self):
        return self.variables

    def get_variables(self, requested_variables, time=None,
                  x=None, y=None, z=None, block=False):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)
        
        nearestTime, tidx, tidx_nxt, indxTime, time_weight, _ = self.nearest_time(time)
        variables = {}
        
        # Determinar nivel de profundidad
        z_level = 0
        if self.z is not None and z is not None:
            valid_z = -np.abs(z)  # Profundidad como valores negativos
            z_idx = np.argmin(np.abs(self.z - valid_z))
            z_level = self.z[z_idx]

        for par in requested_variables:
            # Encontrar variable NEMO correspondiente
            varname = next((k for k, v in self.NEMO_variable_mapping.items() if v == par), None)
            if not varname:
                logging.warning(f"Variable {par} no encontrada en el mapping.")
                continue

            # Seleccionar dataset apropiado
            ds = self.ds
            mask = self.tmask
            if varname == 'uoce' and self.dsU is not None:
                ds = self.dsU
                mask = self.umask
            elif varname == 'voce' and self.dsV is not None:
                ds = self.dsV
                mask = self.vmask
            elif varname == 'wo' and self.dsW is not None:
                ds = self.dsW
            elif varname == 'tmask' and self.ds_mask is not None:
                ds = self.ds_mask

            # Obtener datos 3D (tiempo, lat, lon)
            try:
                if varname in ds:
                    var_data = ds[varname]
                    
                    # Manejar diferentes dimensiones
                    if 'deptht' in var_data.dims and z is not None:
                        data_slice = var_data.isel(time_counter=tidx, deptht=z_idx)
                    elif 'time_counter' in var_data.dims:
                        data_slice = var_data.isel(time_counter=tidx)
                    else:
                        data_slice = var_data
                    
                    # Convertir a numpy array y aplicar máscara
                    raw_data = data_slice.values
                    if mask is not None and raw_data.ndim >= 2:
                        if mask.ndim == 3:  # (time, lat, lon)
                            mask_slice = mask[tidx, :, :]
                        elif mask.ndim == 2:  # (lat, lon)
                            mask_slice = mask
                        # Asegurar que las dimensiones coincidan
                        if raw_data.shape == mask_slice.shape:
                            raw_data = np.where(mask_slice > 0.5, raw_data, np.nan)
                    
                    # Interpolar a las posiciones de las partículas
                    data = scipy.ndimage.map_coordinates(
                        raw_data,
                        np.array([
                            np.interp(y, self.yi, np.arange(len(self.yi))),
                            np.interp(x, self.xi, np.arange(len(self.xi)))
                        ]),
                        order=1,
                        cval=np.nan
                    )
                    
                    # Manejar valores faltantes
                    data = np.ma.masked_invalid(data)
                    data = np.ma.masked_outside(data, -10, 10)  # Rango físico realista

                else:
                    raise ValueError(f"Variable {varname} no encontrada en dataset")
                    
            except Exception as e:
                logging.error(f"Error procesando {par}: {str(e)}")
                data = np.ma.masked_all(len(x))
            
            # Formatear salida para OpenDrift
            data = np.ma.array(np.atleast_2d(data), dtype=np.float32)
            variables[par] = data

        # Añadir coordenadas y tiempo
        variables['z'] = z_level * np.ones(len(x))
        variables['x'] = x
        variables['y'] = y
        variables['time'] = nearestTime

        return variables