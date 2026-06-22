import logging
import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from pyproj import CRS
from opendrift.readers.basereader import BaseReader

class Reader(BaseReader):

    def __init__(self, filename=None, filenameU=None, filenameV=None,
                 filenameW=None, filename_mask=None, name=None,
                 gridfile=None, custom_var_mapping=None, _FillValue=None):

        if filename is None:
            raise ValueError('Necesitas proporcionar filename (grid T)')

        # Mapeo de tus variables
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

        try:
            # Archivos T,U,V,W por separado
            self.Dataset = Dataset(filename, 'r')
            self.Dataset.set_auto_maskandscale(False)

            if filenameU is not None:
                self.DatasetU = Dataset(filenameU, 'r')
                self.DatasetU.set_auto_maskandscale(False)
            else:
                self.DatasetU = None

            if filenameV is not None:
                self.DatasetV = Dataset(filenameV, 'r')
                self.DatasetV.set_auto_maskandscale(False)
            else:
                self.DatasetV = None

            if filenameW is not None:
                self.DatasetW = Dataset(filenameW, 'r')
                self.DatasetW.set_auto_maskandscale(False)
            else:
                self.DatasetW = None

            if filename_mask is not None:
                self.Dataset_mask = Dataset(filename_mask, 'r')
            else:
                self.Dataset_mask = None

        except Exception as e:
            raise ValueError(f"Error abriendo archivos: {e}")

        self._FillValue = _FillValue

        # Coordenadas (preferible desde gridfile o desde Dataset T)
        if gridfile:
            gf = Dataset(gridfile)
            self.lat = gf.variables['nav_lat'][:]
            self.lon = gf.variables['nav_lon'][:]
        else:
            self.lat = self.Dataset.variables.get('nav_lat', None)
            self.lon = self.Dataset.variables.get('nav_lon', None)
            if self.lat is not None:
                self.lat = self.lat[:]
            if self.lon is not None:
                self.lon = self.lon[:]

        # self.projection = CRS.from_epsg(4326)
        self.proj4 = CRS.from_epsg(4326).to_proj4()
        # Tiempo (desde archivo T)
        ocean_time = self.Dataset.variables['time_counter']
        time_units = ocean_time.units
        self.times = num2date(ocean_time[:], time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.time_step = (self.times[1] - self.times[0]) if len(self.times) > 1 else None

        # Profundidad (buscar en Dataset T)
        for dvn in ['deptht', 'depthu', 'depthv', 'depthw']:
            if dvn in self.Dataset.variables:
                varz = self.Dataset.variables[dvn]
                if 'positive' not in varz.ncattrs() or varz.getncattr('positive') == 'up':
                    self.z = varz[:]
                else:
                    self.z = -varz[:]
                break
        else:
            self.z = None

        # Dimensiones espaciales
        self.xmin = 0
        self.xmax = len(self.Dataset.dimensions['x']) - 1
        self.delta_x = 1
        self.ymin = 0
        self.ymax = len(self.Dataset.dimensions['y']) - 1
        self.delta_y = 1

        # Variables disponibles
        self.variables = []
        for var_name in self.NEMO_variable_mapping:
            self.variables.append(self.NEMO_variable_mapping[var_name])

        super().__init__()
        
    @property
    def variables_available(self):
        return self.variables

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        import logging

        # Validación y ajustes iniciales
        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, _, _, indxTime, _, _ = self.nearest_time(time)
        variables = {}

        # Índice vertical
        if self.z is not None and z is not None:
            indices = np.searchsorted(-self.z, [-z.min(), -z.max()])
            indz = np.arange(max(0, indices.min() - 1),
                             min(len(self.z), indices.max() + 1))
            if len(indz) == 1:
                indz = indz[0]
        else:
            indz = 0

        # Convertir coordenadas lon/lat a índices reales de malla
        indx = np.zeros_like(x, dtype=int)
        indy = np.zeros_like(y, dtype=int)
        for i in range(len(x)):
            ii, jj = self.lonlat_to_ij(x[i], y[i])
            indx[i] = ii
            indy[i] = jj

        # Enmascarar partículas fuera del dominio
        indx[outside] = 0
        indy[outside] = 0

        for par in requested_variables:
            varnames = [k for k, v in self.NEMO_variable_mapping.items() if v == par]
            if not varnames:
                logging.warning(f"Variable {par} no encontrada en mapping")
                continue
            varname = varnames[0]

            # Dataset según variable
            if varname == 'uoce' and self.DatasetU is not None:
                var = self.DatasetU.variables[varname]
            elif varname == 'voce' and self.DatasetV is not None:
                var = self.DatasetV.variables[varname]
            elif varname == 'wo' and self.DatasetW is not None:
                var = self.DatasetW.variables[varname]
            elif varname == 'tmask' and self.Dataset_mask is not None:
                var = self.Dataset_mask.variables[varname]
            else:
                var = self.Dataset.variables[varname]

            var.set_auto_maskandscale(False)
            FillValue = getattr(var, '_FillValue', self._FillValue)
            scale = getattr(var, 'scale_factor', 1)
            offset = getattr(var, 'add_offset', 0)

            shape = var.shape
            time_len = shape[0] if len(shape) > 0 else 0
            depth_len = shape[1] if len(shape) > 1 else 0
            lat_len = shape[2] if len(shape) > 2 else 0
            lon_len = shape[3] if len(shape) > 3 else 0

            indxTime_safe = min(indxTime, time_len - 1) if time_len > 0 else 0
            if isinstance(indz, np.ndarray):
                indz_safe = indz[indz < depth_len]
                if indz_safe.size == 0:
                    indz_safe = 0
            else:
                indz_safe = min(indz, depth_len - 1) if depth_len > 0 else 0

            indy_safe = np.clip(indy, 0, lat_len - 1) if lat_len > 0 else indy
            indx_safe = np.clip(indx, 0, lon_len - 1) if lon_len > 0 else indx

            try:
                if var.ndim == 2:
                    data = var[indy_safe, indx_safe]
                elif var.ndim == 3:
                    data = var[indxTime_safe, indy_safe, indx_safe]
                elif var.ndim == 4:
                    data = var[indxTime_safe, indz_safe, indy_safe, indx_safe]
                else:
                    raise Exception(f'Dimensión incorrecta para variable {par}: {var.shape}')
            except Exception as e:
                logging.error(f'Error al leer variable {par}: {e}')
                continue

            data = data * scale + offset

            if FillValue is not None:
                mask = data == FillValue
                data = np.ma.array(data, mask=mask)

            if data.ndim > 1 and not block:
                data = data.diagonal()

            data = np.ma.array(data.copy(), ndmin=2, mask=False)
            data.mask[outside] = True
            data = np.ma.masked_outside(data, -30000, 30000)

            variables[par] = data

        variables['z'] = self.z[indz_safe] if self.z is not None else None
        variables['x'] = x
        variables['y'] = y
        variables['time'] = nearestTime

        return variables

    
    def lonlat_to_ij(self, lon, lat):
        """
        Convierte coordenadas lon/lat a índices i, j en la malla del modelo.
    
        Parámetros:
            lon: longitud (float)
            lat: latitud (float)
    
        Retorna:
            (i, j): índices enteros más cercanos en la grilla
        """
        if self.lon is None or self.lat is None:
            raise ValueError("La grilla nav_lon/nav_lat no fue cargada.")
    
        dist2 = (self.lon - lon) ** 2 + (self.lat - lat) ** 2
        j, i = np.unravel_index(np.argmin(dist2), self.lon.shape)
        return i, j
