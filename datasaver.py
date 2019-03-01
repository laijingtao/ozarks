import numpy as np
from netCDF4 import Dataset

class DataSaver(object):
    def __init__(self, nrows, ncols, dx, dt):
        self.nrows = nrows
        self.ncols = ncols
        self.nt = 0
        self.dx = dx
        self.dt = dt
        self.data = {}
        self.data_units = {}

    def add_field(self, field, unit=None):
        self.data[field] = np.array([])
        self.data_units[field] = unit

    def add_value(self, field, z):
        if field not in self.data:
            raise KeyError(field)
        if len(self.data[field]) == 0:
            self.data[field] = np.zeros((1, self.nrows, self.ncols))
            self.data[field][0] = np.array(z).reshape((self.nrows, self.ncols))
        else:
            self.data[field] = np.append(self.data[field],
                                         [np.array(z).reshape((self.nrows, self.ncols))], axis=0)
        if len(self.data[field]) > self.nt:
            self.nt = len(self.data[field])

    def write_to_nc(self, file_name=None):
        if file_name is None:
            raise ValueError("File name")

        ncdata = Dataset(file_name, 'w')

        dim_value = {'x': np.arange(self.ncols)*self.dx,
                     'y': np.arange(self.nrows)*self.dx,
                     'time': np.arange(self.nt)*self.dt}
        dim_len = {'x': self.ncols, 'y': self.nrows, 'time': self.nt}
        for dim_name in ['x', 'y', 'time']:
            try:
                ncdata.createDimension(dim_name, dim_len[dim_name])
                ncdata.createVariable(dim_name, np.float64, (dim_name,))
                if dim_len[dim_name] is not None:
                    ncdata.variables[dim_name][:] = dim_value[dim_name]
            except:
                continue
        ncdata.variables['time'].units = 'model years'
        #ncdata.variables['time'].calendar = 'standard'

        for var_name in self.data:
            try:
                var = ncdata.variables[var_name]
            except:
                var = ncdata.createVariable(var_name, np.float64, ('time', 'y', 'x',))
            #print(var.shape)
            var[:] = self.data[var_name]
            if self.data_units[var_name] is not None:
                var.units = self.data_units[var_name]

        ncdata.close()
