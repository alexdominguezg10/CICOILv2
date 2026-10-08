"""
OpenDrift reader for CROCO native output (croco_avg.nc / croco_his.nc), no preprocessing.

Subclass of OpenDrift 1.14.9's reader_ROMS_native, which already handles CROCO naming (Cs_rho, time,
angle, lon_rho, sigma -> z on the fly). Two things are added for CROCO files:

 1. mask_u / mask_v are derived from mask_rho (CROCO files only carry mask_rho): a u- or v-point is wet only
    if both neighbouring rho-points are wet.
 2. destagger=True (default): the stock reader skips de-staggering and treats u (xi_u) and v (eta_v) values
    as if they sat on the rho nodes, i.e. a half-cell (~2.3 km on the 1/24 deg grid) position error.
    Here u and v are averaged onto the rho points lazily (dask) right after opening, with the end points
    copied -- the same operation as in prep_forcing.py -- so the reader then works on a plain rho grid.
    destagger=False reproduces the stock behaviour.
"""
import logging
import numpy as np, xarray as xr
import dask.array as da
import opendrift
from opendrift.readers.reader_ROMS_native import Reader as ROMSReader

logger = logging.getLogger(__name__)
TESTED_OPENDRIFT = '1.14.9'


def _pad_avg(arr, ax):
    n = arr.shape[ax]

    def s(a, b):
        idx = [slice(None)] * arr.ndim
        idx[ax] = slice(a, b)
        return arr[tuple(idx)]
    return da.concatenate([s(0, 1), 0.5 * (s(0, n - 1) + s(1, n)), s(n - 1, n)], axis=ax)


class Reader(ROMSReader):
    def __init__(self, *args, destagger=True, **kwargs):
        if opendrift.__version__ != TESTED_OPENDRIFT:
            logger.warning(f'reader_croco_native was written and tested with OpenDrift {TESTED_OPENDRIFT} '
                           f'(found {opendrift.__version__}); it relies on private attributes of reader_ROMS_native.')
        super().__init__(*args, **kwargs)
        self.destaggered = False
        if destagger:
            self._destagger()

    def _destagger(self):
        ds = self.Dataset
        for name, stag in (('u', 'xi_u'), ('v', 'eta_v')):
            if name not in ds.variables:
                continue
            v = ds[name]
            ax = v.dims.index(stag)
            data = v.data if hasattr(v.data, 'chunks') else da.from_array(np.asarray(v.data), chunks='auto')
            data = da.nan_to_num(data)
            new = _pad_avg(data, ax)
            dims = tuple({'eta_u': 'eta_rho', 'xi_u': 'xi_rho', 'eta_v': 'eta_rho', 'xi_v': 'xi_rho'}.get(d, d) for d in v.dims)
            ds[name] = xr.DataArray(new, dims=dims, attrs=v.attrs)
        self.destaggered = True

    @property
    def mask_u(self):
        if self._mask_u is None:
            m = np.asarray(self.mask_rho)
            self._mask_u = xr.DataArray(m[:, 1:] * m[:, :-1], dims=('eta_u', 'xi_u'))
        return self._mask_u

    @property
    def mask_v(self):
        if self._mask_v is None:
            m = np.asarray(self.mask_rho)
            self._mask_v = xr.DataArray(m[1:, :] * m[:-1, :], dims=('eta_v', 'xi_v'))
        return self._mask_v
