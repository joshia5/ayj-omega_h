# In case omega_h is not in the default path, append the path where the PyOmega_h shared library is located
import sys
sys.path.append(
  "/lore/joshia5/develop/mfem_omega/build-omegah-python-rhel7/src")
import PyOmega_h as omega_h

#import ctypes
#from ctypes import cdll, c_char_p
#libc = cdll.LoadLibrary("libc.so.6")
#puts = libc.puts
#puts('abc')
#puts.argtypes = [c_char_p]
#puts(b'abc')
infile = "/lore/joshia5/Meshes/curved/inclusion_3p_sizes_2.osh"
#b_infile = infile.encode('utf-8')
#puts(b_infile)
#puts(infile.encode())
#path = omega_h.path_c(b_infile)

# Build mesh with omega_h & initialize metrics
comm_osh = omega_h.world()
#mesh_osh = omega_h.build_box(comm_osh, omega_h.Family.SIMPLEX, 1.0, 1.0, 1.0, 8, 8, 8)
#omega_h.vtk_write_parallel('test_adapted_' + str(0), mesh_osh)
#mesh_osh.balance()
#mesh_osh.set_parting(omega_h.GHOSTED, 0)
#omega_h.add_implied_metric_tag(mesh_osh)
#mesh_osh.set_parting(omega_h.ELEM_BASED, 0)
#binary_read_file = omega_h.binary_read_file
#mesh_osh = omega_h.new_empty_mesh()
#binary_read_file.argtypes = [c_char_p, omega_h.Comm, bool]
#mesh_osh = omega_h.gmsh_read_file(ctypes.c_char_p(b_infile), comm_osh);
mesh_osh = omega_h.binary_read_file(infile, comm_osh);
#mesh_osh = binary_read_file(ctypes.c_char_p(b_infile), comm_osh);
#mesh_osh = omega_h.binary_read_file(ctypes.c_char_p(b_infile), comm_osh);
#mesh_osh = omega_h.binary_read_file(infile.encode()._cpp_object, comm_osh);
#mesh_osh = omega_h.binary_read_file('/lore/joshia5/Meshes/curved/inclusion_3p_sizes.osh', comm_osh);
#mesh_osh = omega_h.binary_read_file(ctypes.c_char_p("/lore/joshia5/Meshes/curved/inclusion_3p_sizes.osh"), comm_osh);

maxiter = 2
i = 0

while(i < maxiter):
 
    #mesh_osh.set_parting(omega_h.GHOSTED, 0);
    #metric_input = omega_h.MetricInput()
    #source = omega_h.MetricSource(omega_h.VARIATION, 2e-3, "u")
    #metric_input.add_source(source)
    #metric_input.should_limit_lengths = True
    #metric_input.max_length = 0.1
    #metric_input.should_limit_gradation = True
    print("mesh nents v e f r\n", mesh_osh.nents(0),mesh_osh.nents(1),
        mesh_osh.nents(2), mesh_osh.nents(3))
    #omega_h.add_implied_metric_tag(mesh_osh)
    #omega_h.generate_target_metric_tag(mesh_osh, metric_input) 

    opts = omega_h.AdaptOpts(mesh_osh)
    opts.should_swap = 0
    opts.should_coarsen = 0
    opts.should_coarsen_slivers = 0
    opts.should_filter_invalids = 0
    opts.check_crv_qual = 0
    opts.verbosity = omega_h.WRITE_FILE
    #mesh_osh.add_tag<omega_h.Real>(0, "metric", 1);
    #mesh_osh.set_tag(0, "metric", omega_h.Reals(mesh_osh.nverts(),
      #omega_h.metric_eigenvalue_from_length(100)));

    # Adapt mesh
    while(omega_h.approach_metric(mesh_osh, opts) and
      mesh_osh.nents(3) < 10000):
        omega_h.adapt(mesh_osh, opts)

    #omega_h.vtk_write_parallel('vtk_test_adapted_' + str(i), mesh_osh)
    #omega_h.binary_write_file('test_adapted_' + str(i), mesh_osh)

    i+=1

