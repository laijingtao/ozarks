from __future__ import print_function
from base import *

def run_model(params_file=None):

    if params_file is None:
        params_file = 'params_test.txt'
    params = ParamParser(params_file)

    runtime = int(params.read('runtime', 'float')) # years
    dt = params.read('dt', 'int') # years
    nt = int(runtime/dt)
    Lx = params.read('Lx', 'float')
    Ly = params.read('Ly', 'float')
    dx = params.read('dx', 'float') # meter
    nrows = int(Ly/dx)
    ncols = int(Lx/dx)
    #nrows = params.read('nrows', 'int')
    #ncols = params.read('ncols', 'int')
    num_of_nodes = nrows*ncols
    #uplift_rate = 0.0015 # m/year
    k_sp_orogen = params.read('k_sp_orogen', 'float')
    k_sp_foreland = params.read('k_sp_foreland', 'float')
    elastic_thickness = params.read('elastic_thickness', 'float')
    save_dt = params.read('save_dt', 'int')
    num_plots = params.read('num_plots', 'int')
    if num_plots == 0:
        plot_interval = runtime + dt
    else:
        plot_interval = int(runtime/num_plots)
    outdir = params.read('outdir', 'str')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    #if not os.path.isdir(os.path.join(outdir, 'tmp')):
    #    os.mkdir(os.path.join(outdir, 'tmp'))
    outfile = params.read('outfile', 'str')
    outfile_base, _ = os.path.splitext(outfile)

    try:
        initial_file = params.read('initial_file', 'str')
    except MissingKeyError:
        initial_file = 'initial_topg.nc'
    initial_mg = read_netcdf(initial_file)

    mg = RasterModelGrid((nrows, ncols), dx)

    orogen_width = params.read('orogen_width', 'int')
    foreland = np.intersect1d(mg.core_nodes, np.where(mg.node_y > orogen_width)[0])  # foreland basin
    orogen = np.intersect1d(mg.core_nodes, np.where(mg.node_y <= orogen_width)[0])  # mountain, load
    #river_outlet = np.intersect1d(np.where(mg.node_y == orogen_width+dx)[0],
    #                              np.where(mg.node_x == mg.node_x.max())[0])  # base river
    river_outlet = np.intersect1d(np.where(mg.node_y == orogen_width+dx)[0],
                                  np.where(np.logical_or(mg.node_x == mg.node_x.max(),
                                                         mg.node_x == mg.node_x.min()))[0])  # base river
    river = np.where(mg.node_y == orogen_width+dx)[0]

    mg.add_zeros('node', 'topographic__elevation', units='m')
    z = mg.at_node['topographic__elevation']
    z[:] = initial_mg.at_node['topographic__elevation'][:]
    #z[foreland] += np.random.rand(len(z[foreland])) # may not need this if initial has a drainage network

    mg.add_zeros('node', 'stream_power_k_field')
    k_field = mg.at_node['stream_power_k_field']
    k_field[foreland] = k_sp_foreland
    k_field[orogen] = k_sp_orogen
    k_field[river] = k_sp_orogen

    diff_field = k_field/0.002

    #set up grid's boundary conditions (right, top, left, bottom) is inactive
    #mg.set_closed_boundaries_at_grid_edges(True, False, True, True)
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)
    #mg.status_at_node[river_outlet] = FIXED_VALUE_BOUNDARY
    #mg.status_at_node[river] = FIXED_VALUE_BOUNDARY

    fr = FlowRouter(mg)
    sp = FastscapeEroder(mg, K_sp=k_field)
    lin_diffuse = LinearDiffuser(mg, linear_diffusivity=diff_field)

    beam_length = params.read('beam_length', 'float')
    gflex_mg = RasterModelGrid(5, (int(beam_length/dx)), dx)
    gflex_mg.add_zeros('node', 'topographic__elevation', units='m')
    gflex_mg.add_zeros('node', 'surface_load__stress')
    load = gflex_mg.at_node['surface_load__stress']

    # it seems BC_S and BC_N are reversed in landlab's gflex
    gf = gFlex(gflex_mg, elastic_thickness=elastic_thickness,
               BC_S='Mirror', BC_N='Mirror', BC_W='Mirror', BC_E='0Displacement0Slope')

    def update_flexure(mg, gflex_mg, gf, topg_change):
        topg_change = mg.node_vector_to_raster(topg_change)
        mean_topg_change = np.array([np.mean(topg_change[k, :]) for k in range(nrows)])
        load[np.where(gflex_mg.node_x < nrows*dx)] = np.array([-9.8*2700*mean_topg_change for k in range(5)]).reshape(5*nrows)
        gf.run_one_step()
        displacement = gflex_mg.at_node['lithosphere_surface__elevation_increment']
        displacement_1d = displacement[np.where(gflex_mg.node_y == 2*dx)][0:nrows]
        displacement_2d = np.zeros((nrows, ncols))
        for k in range(ncols):
            displacement_2d[:, k] = displacement_1d
        displacement_2d = displacement_2d.reshape(nrows*ncols)
        mg.at_node['topographic__elevation'][mg.core_nodes] += displacement_2d[mg.core_nodes]
        gflex_mg.at_node['topographic__elevation'][:] = 0
        gflex_mg.at_node['lithosphere_surface__elevation_increment'][:] = 0

    mg.add_zeros('node', 'erosion__accum', units='m')
    mg.add_zeros('node', 'uplift__accum', units='m')

    '''
    write_netcdf(os.path.join(outdir, 'tmp', '{}_t=0.nc'.format(outfile_base)),
                 mg, format='NETCDF3_64BIT',
                 names=['topographic__elevation', 'erosion__accum', 'uplift__accum'])
    '''
    write_raster_netcdf(os.path.join(outdir, outfile), mg, time=0, append=False,
                        names=['topographic__elevation', 'erosion__accum', 'uplift__accum'])

    previous_topg = np.zeros(len(z))
    previous_topg[:] = z[:]
    topg = mg.node_vector_to_raster(z)
    previous_mean_topg = np.array([np.mean(topg[k, :]) for k in range(nrows)])
    topg_change_erosion = np.zeros(len(z))
    topg_change_uplift = np.zeros(len(z))

    plot_idx = 1
    save_interval = 0
    plt.close('all')
    for i in range(nt):
        fr.run_one_step()
        previous_topg[:] = z[:]
        sp.run_one_step(dt)
        topg_change_erosion = previous_topg - z
        mg.at_node['erosion__accum'] += topg_change_erosion
        previous_topg[:] = z[:]

        # uplift rate is the highest on the bottom boundary, decreases to 0 linearly on the top boundary
        #z[mg.core_nodes] += uplift_rate * dt * (mg.node_y.max() - mg.node_y[mg.core_nodes]) / mg.node_y.max()

        update_flexure(mg, gflex_mg, gf, topg_change_erosion)
        topg_change_uplift = z - previous_topg
        mg.at_node['uplift__accum'] += topg_change_uplift
        previous_topg[:] = z[:]

        lin_diffuse.run_one_step(dt)

        save_interval += dt
        if save_interval == save_dt:
            '''
            write_netcdf(os.path.join(outdir, 'tmp', '{}_t={}.nc'.format(outfile_base, int((i+1)*dt))),
                         mg, format='NETCDF3_64BIT',
                         names=['topographic__elevation', 'erosion__accum', 'uplift__accum'])
            '''
            write_raster_netcdf(os.path.join(outdir, outfile), mg, time=(i+1)*dt, append=True,
                                names=['topographic__elevation', 'erosion__accum', 'uplift__accum'])
            save_interval = 0

        if int(dt*(i+1)) == plot_interval*plot_idx:
            fig = plt.figure(plot_idx, figsize=(8, 9))
            imshow_grid(mg, z, cmap='terrain', grid_units=('m', 'm'),
                        plot_name='t = {} kyr'.format(int(dt*(i+1)/1000.0)), var_name='Elevation (m)')

            #plt.figure(plot_idx + num_plots)
            ax = fig.add_axes([1.1, 0.2, 0.6, 0.6])
            y = mg.node_vector_to_raster(mg.node_y)[1:nrows, 0]
            topg = mg.node_vector_to_raster(z)
            plt_mean_topg = [np.mean(topg[k, :]) for k in range(1, nrows)]
            ax.plot(y, plt_mean_topg)
            ax.set_xlabel('Y (m)')
            ax.set_ylabel('Elevation (m)')
            ax.set_title('Mean elevation (t = {} kyr)'.format(int(dt*(i+1)/1000.0)))
            plot_idx += 1

        #print('Running... [{}%]'.format(int((i+1)*100.0/nt)), end='\r')

    #post_process(params_file=params_file)

def post_process(params_file=None):
    if params_file is None:
        params_file = 'params_test.txt'
    params = ParamParser(params_file)

    runtime = int(params.read('runtime', 'float')) # years
    nrows = params.read('nrows', 'int')
    ncols = params.read('ncols', 'int')
    dx = params.read('dx', 'float') # meter
    save_dt = params.read('save_dt', 'int')
    nt = int(runtime/save_dt)
    outdir = params.read('outdir', 'str')
    outfile = params.read('outfile', 'str')
    outfile_base, _ = os.path.splitext(outfile)

    datasaver = DataSaver(nrows, ncols, dx, save_dt)
    datasaver.add_field('topographic__elevation', 'm')
    datasaver.add_field('erosion__accum', 'm')
    datasaver.add_field('uplift__accum', 'm')

    mg = read_netcdf(os.path.join(outdir, 'tmp', '{}_t=0.nc'.format(outfile_base)))
    datasaver.add_value('topographic__elevation', mg.at_node['topographic__elevation'])
    datasaver.add_value('erosion__accum', mg.at_node['erosion__accum'])
    datasaver.add_value('uplift__accum', mg.at_node['uplift__accum'])

    for i in range(nt):
        mg = read_netcdf(os.path.join(outdir, 'tmp', '{}_t={}.nc'.format(outfile_base, int((i+1)*save_dt))))
        datasaver.add_value('topographic__elevation', mg.at_node['topographic__elevation'])
        datasaver.add_value('erosion__accum', mg.at_node['erosion__accum'])
        datasaver.add_value('uplift__accum', mg.at_node['uplift__accum'])

    datasaver.write_to_nc(os.path.join(outdir, outfile))


if __name__ == '__main__':
    run_model()
