###
# Create initial topography
###

from __future__ import print_function
from base import *

params_file = 'params.txt'
params = ParamParser(params_file)

runtime = 5e6 # years
dt = 10000 # years
nt = int(runtime/dt)
nrows = params.read('nrows', 'int')
ncols = params.read('ncols', 'int')
num_of_nodes = nrows*ncols
dx = params.read('dx', 'float') # meter
uplift_rate = 0.0015 # m/year
k_sp_orogen = params.read('k_sp_orogen', 'float')

mg = RasterModelGrid((nrows, ncols), dx)
mg.add_zeros('node', 'topographic__elevation', units='m')
z = mg.at_node['topographic__elevation']
z += np.random.rand(len(z))

orogen_width = params.read('orogen_width', 'int')
foreland = np.intersect1d(mg.core_nodes, np.where(mg.node_y > orogen_width)[0])  # foreland basin
orogen = np.intersect1d(mg.core_nodes, np.where(mg.node_y <= orogen_width)[0])  # mountain, load
river_outlet = np.intersect1d(np.where(mg.node_y == orogen_width)[0],
                              np.where(mg.node_x == mg.node_x.max())[0])  # base river

#set up grid's boundary conditions (right, top, left, bottom) is inactive
mg.set_closed_boundaries_at_grid_edges(True, False, True, True)
mg.status_at_node[river_outlet] = FIXED_VALUE_BOUNDARY

fr = FlowRouter(mg)
sp = FastscapeEroder(mg, K_sp=k_sp_orogen)

for i in range(nt):
    mg.at_node['topographic__elevation'][orogen] += uplift_rate*dt
    fr.run_one_step()
    sp.run_one_step(dt)

    print('Running... [{}%]'.format(int((i+1)*100.0/nt)), end='\r')

write_netcdf('initial_topg.nc', mg, format='NETCDF3_64BIT', names='topographic__elevation')
