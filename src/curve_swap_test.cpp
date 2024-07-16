#include <Omega_h_build.hpp>
#include <Omega_h_swap.hpp>
#include <Omega_h_library.hpp>
//#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
//#include <Omega_h_for.hpp>
//#include <Omega_h_curve_validity_3d.hpp>

using namespace Omega_h;

/*
void test_disc_collapse(Library *lib) {
  auto comm = lib->world();

  //auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_2tri_order2.sms",
    //                    "/users/joshia5/Models/curved/disk_semi_geomsim.smd", comm);
  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
                            "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", comm);
                            
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    //"/users/joshia5/Meshes/curved/disc2_cubic_wireframe.vtu";
    "/users/joshia5/Meshes/curved/disc100_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  //vtuPath = "/users/joshia5/Meshes/curved/disc2_cubic_curveVtk.vtu";
  vtuPath = "/users/joshia5/Meshes/curved/disc100_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  auto coords = mesh.coords();
  auto ev2v = mesh.get_adj(1, 0).ab2b;
  auto fe2e = mesh.get_adj(2, 1).ab2b;
  auto v_class_dim = mesh.get_array<I8>(0, "class_dim");
  auto v_class_id = mesh.get_array<LO>(0, "class_id");
  auto e_class_dim = mesh.get_array<I8>(1, "class_dim");
  auto e_class_id = mesh.get_array<LO>(1, "class_id");
  printf("coords %d, ev2v %d, fe2e %d\n", coords.size(), ev2v.size(), fe2e.size());
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(0.5)));
  while ((coarsen_by_size(&mesh, opts)) && (mesh.nelems() > 10));
  mesh.ask_qualities();
  coords = mesh.coords();
  ev2v = mesh.get_adj(1, 0).ab2b;
  fe2e = mesh.get_adj(2, 1).ab2b;
  return;
}


void test_collapse_boxCircle(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "../omega_h/meshes/box_circleCut_4k.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  vtk::FullWriter writer;

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "../omega_h/meshes/box_circleCut_4k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "../omega_h/meshes/box_circleCut_4k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());
  I8 max_adapt_itr = 4;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    coarsen_by_size(&mesh, opts);
    //coarsen_slivers(&mesh, opts);
  }
  mesh.ask_qualities();
  writer = vtk::FullWriter(
      "../omega_h/meshes/boxCircle_aft.vtk", &mesh);
  writer.write();
  return;
}

void test_antenna(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/antenna_6k.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/antenna_6k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/antenna_6k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  coarsen_by_size(&mesh, opts);
  mesh.ask_qualities();
  return;
}

void test_collapse_cubicSlab(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/cubic_slab-case1.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/cubic_slab_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/cubic_slab_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
 
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.00001;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(10000)));
  coarsen_by_size(&mesh, opts);
  mesh.ask_qualities();
  return;
}

void test_sphere(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/lore/joshia5/Meshes/curved/sphere_8.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/sphere_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/sphere_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
 
  //binary::write("/lore/joshia5/Meshes/curved/sphere_8_invalid.osh", &mesh);
  
  //auto opts = AdaptOpts(&mesh);
  //opts.min_quality_allowed = 0.00001;
  //mesh.add_tag<Real>(VERT, "metric", 1);
  //mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
  //      metric_eigenvalue_from_length(10000)));
  //coarsen_by_size(&mesh, opts);
  //mesh.ask_qualities();
  
  return;
}
*/

void test_swap_kova(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("../omega_h/meshes/KovaGeomSim-quadratic_123tet.osh", comm);
  printf("initial ntets %d\n", mesh.nregions());
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  vtk::FullWriter writer;
  writer = vtk::FullWriter("/lore/joshia5/Meshes/curved/kovaCoarsen_bef.vtk", &mesh);
  writer.write();
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);

  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
 
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 5);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/Kova_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 5);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.9;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
        //metric_eigenvalue_from_length(1)));
  fprintf(stderr, "start swapping\n");
  for (LO adapt_itr = 0; adapt_itr < 3; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }
  fprintf(stderr, "finish swapping\n");
  writer = vtk::FullWriter("/lore/joshia5/Meshes/curved/kovaSwap_aft.vtk", &mesh);
  writer.write();
  fprintf(stderr, "finished kova case\n");
  
  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_swap_kova(&lib);
  //test_disc_collapse(&lib);
  //test_disc_validity(&lib);
  //test_collapse_boxCircle(&lib);
  //test_collapse_cubicSlab(&lib);
  //test_antenna(&lib);
  //test_sphere(&lib);

  return 0;
}
