# In case omega_h is not in the default path, append the path where the PyOmega_h shared library is located
import sys
sys.path.append("/lore/joshia5/develop/mfem_omega/build-omegah-python-rhel7/src")
import PyOmega_h as omega_h;

# Build mesh with omega_h & initialize metrics
comm_osh = omega_h.world()
mesh_osh = omega_h.build_box(comm_osh, omega_h.Family.SIMPLEX, 1.0, 1.0, 1.0, 8, 8, 8)
mesh_osh.balance()
mesh_osh.set_parting(omega_h.GHOSTED, 0)
omega_h.add_implied_metric_tag(mesh_osh)
mesh_osh.set_parting(omega_h.ELEM_BASED, 0)

maxiter = 2
i = 0

while(i < maxiter):
 
    #mesh_osh.set_parting(omega_h.GHOSTED, 0);
    metric_input = omega_h.MetricInput()
    #source = omega_h.MetricSource(omega_h.VARIATION, 2e-3, "u")
    #metric_input.add_source(source)
    metric_input.should_limit_lengths = True
    #metric_input.max_length = 0.1
    #metric_input.should_limit_gradation = True
    print("ok26\n")
    omega_h.add_implied_metric_tag(mesh_osh)
    print("ok27\n")
    omega_h.generate_target_metric_tag(mesh_osh, metric_input) 
    print("ok28\n")

    opts = omega_h.AdaptOpts(mesh_osh)
    opts.verbosity = omega_h.EXTRA_STATS
    #mesh_osh.add_tag<omega_h.Real>(0, "metric", 1);
    #mesh_osh.set_tag(0, "metric", omega_h.Reals(mesh_osh.nverts(),
      #omega_h.metric_eigenvalue_from_length(100)));

    # Adapt mesh
    while(omega_h.approach_metric(mesh_osh, opts)):
        omega_h.adapt(mesh_osh, opts)

    omega_h.vtk_write_parallel('test_adapted_' + str(i), mesh_osh)
    
    i+=1

