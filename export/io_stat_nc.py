import sys
from datetime import datetime, timedelta
import string
import logging

import numpy as np
from netCDF4 import Dataset, num2date, date2num

# Module with functions to export/import trajectory data to/from netCDF file
# Strives to be compliant with netCDF CF-convention on trajectories
# http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#idp8377728
# https://geo-ide.noaa.gov/wiki/index.php?title=NODC_NetCDF_Trajectory_Template

def init(self, model_template, filename, prop_metadata):

    self.outfile_name = filename
    self.outfile = Dataset(filename, 'w')

    self.outfile.Conventions = 'CF-1.6'
    self.outfile.standard_name_vocabulary = 'CF-1.6'
    self.outfile.history = 'Created ' + str(datetime.now())
    self.outfile.source = 'Output from simulation with OpenDrift'
    self.outfile.model_url = 'https://github.com/OpenDrift/opendrift'
    self.outfile.opendrift_class = model_template.__class__.__name__
    self.outfile.opendrift_module = model_template.__class__.__module__
    self.outfile.readers = str(model_template.readers.keys())
    self.outfile.time_step_calculation = str(model_template.time_step)
    self.outfile.time_step_output = str(model_template.time_step_output)

    # Write config settings
    for key in model_template._config:
        value = model_template.get_config(key)
        if isinstance(value, (bool, type(None))):
            value = str(value)
        self.outfile.setncattr('config_' + key, value)

    # Create netCDF dimensions
    for prop in prop_metadata:
        if prop_metadata[prop]['type']=='dimension':
            self.outfile.createDimension(prop, prop_metadata[prop]['length'])
            var = self.outfile.createVariable(prop, prop_metadata[prop]['dtype'], prop_metadata[prop]['axes'])
            for subprop in prop_metadata[prop].items():
                if subprop[0] not in ['dtype', 'type', 'axes']:
                    var.setncattr(subprop[0], subprop[1])

    # Create netCDF properties fields
    for prop in prop_metadata:
        if prop_metadata[prop]['type']=='property':
            var = self.outfile.createVariable(prop, prop_metadata[prop]['dtype'], prop_metadata[prop]['axes'])
            for subprop in prop_metadata[prop].items():
                if subprop[0] not in ['dtype', 'type', 'axes']:
                    var.setncattr(subprop[0], subprop[1])
    self.sim_idx = 0
    self.prop_metadata = prop_metadata

# Populate the Dataset with data per simulation
def write_stats_to_nc(self, properties_data, prop_metadata):
    self.yearly = properties_data
    for prop in prop_metadata:
        self.prop = prop
        self.prop_data = properties_data[prop]
        self.axes_no = len(prop_metadata[prop]['axes'])
        axes_names = prop_metadata[prop]['axes']
        self.axes_size = np.zeros(self.axes_no)
        for idx in range(self.axes_no):
            self.axes_size[idx] = self.prop_metadata[axes_names[idx]]['length']
        self.axes_idx = np.zeros(self.axes_no, dtype= np.int32)

        if 'sim_date' in prop_metadata[prop].items(): populate_property(self, self.year)
        else: populate_property(self, self.year)
    self.outfile.sync()  # Flush from memory to disk

def populate_property(self, year=None):
    if self.axes_no == 1:
        if isinstance(self.prop_data, np.ndarray): self.outfile.variables[self.prop][:] = self.prop_data[:]
        else: self.outfile.variables[self.prop][year] = self.prop_data
    elif self.axes_no == 2:
        while self.axes_idx[0] < self.axes_size[0]:
            self.outfile.variables[self.prop][self.axes_idx[0], :] = self.prop_data[self.axes_idx[0],:]
            # print('Time step: ', self.axes_idx[1], '  Values: ', self.prop_data[:, self.axes_idx[1]])
            self.axes_idx[0] +=1
    elif self.axes_no == 3:
        while self.axes_idx[0] < self.axes_size[0]:
            self.outfile.variables[self.prop][self.year, self.axes_idx[0],
                                   self.axes_idx[1], :] = self.prop_data[self.axes_idx[0], self.axes_idx[1], :]
            self.axes_idx[1] +=1
            if self.axes_idx[1] == self.axes_size[1]:
                self.axes_idx[0] +=1
                self.axes_idx[1] = 0
    elif self.axes_no == 4:
        while self.axes_idx[0] < self.axes_size[0]:
            self.outfile.variables[self.prop][self.axes_idx[0], self.axes_idx[1], self.axes_idx[2], :] = self.prop_data[self.axes_idx[0], self.axes_idx[1], self.axes_idx[2], :]
            self.axes_idx[2] += 1
            if self.axes_idx[2] == self.axes_size[2]:
                self.axes_idx[1] += 1
                self.axes_idx[2] = 0
            if self.axes_idx[1] == self.axes_size[1]:
                self.axes_idx[0] += 1
                self.axes_idx[1] = 0

def close(self):
    self.outfile.close()  # Finally close file
