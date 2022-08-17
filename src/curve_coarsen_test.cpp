#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_curve_validity_3d.hpp>

using namespace Omega_h;

void test_tri_validity(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tri.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  auto dim = mesh.dim();
  auto coords = mesh.coords();
  auto edge_nCtrlPts = mesh.n_internal_ctrlPts(1);
  auto ev2v = mesh.ask_down(1, 0).ab2b;

  mesh.add_tag<Real>(0, "bezier_pts", dim, coords);
  mesh.set_tag_for_ctrlPts(1, Reals({1.0/3.0,0.0, 2.0/3.0,0.0,
                                     2.0/3.0,1.0/3.0, 1.0/3.0,2.0/3.0,
                                     0.0,2.0/3.0, 0.0,1.0/3.0}));
  mesh.set_tag_for_ctrlPts(2, Reals({1.0/3.0, 1.0/3.0}));
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/users/joshia5/Meshes/curved/1tri_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/users/joshia5/Meshes/curved/1tri_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris = checkValidity(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  return;
}

/*
void test_disc_validity(Library *lib) {
  auto comm = lib->world();

  
  //auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_2tri_order2.sms",
    //                      "/users/joshia5/Models/curved/disk_semi_geomsim.smd", comm);
  
  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
                            "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", comm);

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/users/joshia5/Meshes/curved/disc2_cubic_wireframe.vtu";
    //"/users/joshia5/Meshes/curved/disc100_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/users/joshia5/Meshes/curved/disc2_cubic_curveVtk.vtu";
  //vtuPath = "/users/joshia5/Meshes/curved/disc100_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  auto valid_tris = checkValidity(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  return;
}

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
*/
void test_linear_tet_validity(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tet.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(1);
  auto vertCtrlPts = HostRead<Real>(Reals({0.0,0.0,0.0,
                                           1.0,0.0,0.0,
                                           0.0,1.0,0.0,
                                           0.0,0.0,1.0}));
  Few<Real, 60> tet_pts;//ntet_pts*dim=20*3
  for (LO j = 0; j < vertCtrlPts.size(); ++j) {
    tet_pts[j] = vertCtrlPts[j];
  }
  /*
  Few<Real, 84> nodes_det = getTetJacDetNodes<84>(1, tet_pts);

  auto is_invalid = checkMinJacDet_3d(nodes_det, 3);
  printf("linear tet is invalid %d \n", is_invalid);
  */
}

void test_quadratic_tet_validity(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tet.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(2);
  auto vertCtrlPts = HostRead<Real>(Reals({0.0,0.0,0.0,
                                           1.0,0.0,0.0,
                                           0.0,1.0,0.0,
                                           0.0,0.0,1.0}));
  auto edgeCtrlPts = HostRead<Real>(Reals({
                                     1.0/2.0, 0.0, 0.0,
                                     //1
                                     1.0/2.0, 1.0/2.0, 0.0,
                                     //3
                                     0.0, 1.0/2.0, 0.0,
                                     //0
                                     0.0, 0.0, 1.0/2.0,
                                     //2
                                     1.0/2.0, 0.0, 1.0/2.0,
                                     //5
                                     0.0, 1.0/2.0, 1.0/2.0
                                     //4
                                     }));
  Few<Real, 60> tet_pts;//ntet_pts*dim=20*3
  for (LO j = 0; j < vertCtrlPts.size(); ++j) {
    tet_pts[j] = vertCtrlPts[j];
  }
  for (LO j = 0; j < edgeCtrlPts.size(); ++j) {
    tet_pts[12 + j] = edgeCtrlPts[j];
  }
  /*
  Few<Real, 84> nodes_det = getTetJacDetNodes<84>(2, tet_pts);

  auto is_invalid = checkMinJacDet_3d(nodes_det, 3);
  printf("quadratic tet is invalid %d \n", is_invalid);
  */
  return;
}

void test_cubic_tet_validity(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tet.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);

  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  auto edge_nCtrlPts = mesh.n_internal_ctrlPts(1);
  auto coords = mesh.coords();
  auto ev2v = mesh.ask_down(1, 0).ab2b;
  auto fe2e = mesh.ask_down(2, 1).ab2b;
  auto rv2v = mesh.ask_down(3, 0).ab2b;
  auto re2e = mesh.ask_down(3, 1).ab2b;
  auto rf2f = mesh.ask_down(3, 2).ab2b;

  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  mesh.set_tag_for_ctrlPts(1, Reals({
                                     0.0, 1.0/3.0, 0.0,
                                     0.0, 2.0/3.0, 0.0,

                                     1.0/3.0, 0.0, 0.0,
                                     2.0/3.0, 0.0, 0.0,

                                     0.0, 0.0, 1.0/3.0,
                                     0.0, 0.0, 2.0/3.0,

                                     1.0/3.0, 2.0/3.0, 0.0,
                                     2.0/3.0, 1.0/3.0, 0.0,

                                     0.0, 1.0/3.0, 2.0/3.0,
                                     0.0, 2.0/3.0, 1.0/3.0,

                                     2.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 0.0, 2.0/3.0
                                     }));
  mesh.set_tag_for_ctrlPts(2, Reals({
                                     1.0/3.0, 1.0/3.0, 0.0,
                                     1.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 1.0/3.0, 1.0/3.0,
                                     0.0, 1.0/3.0, 1.0/3.0
                                     }));
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/tet_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/tet_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  checkValidity_3d(&mesh, LOs({0}));
  
}

void test_boxCircle_validity(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/lore/joshia5/Meshes/curved/box_circleCut-100k.osh", comm);
  if (!mesh.has_tag(0, "bezier_pts")) 
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);

  elevate_curve_order_2to3(&mesh);

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/boxCircle_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/boxCircle_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  checkValidity_3d(&mesh, LOs(mesh.nregions(), 0));
  return;
}

void test_Kova_validity(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet.osh", comm);
  if (!mesh.has_tag(0, "bezier_pts")) 
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);

  elevate_curve_order_2to3(&mesh);

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/Kova_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  checkValidity_3d(&mesh, LOs(mesh.nregions(), 0));
  return;
}

void test_collapse_kova(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read("/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet.osh", comm);
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
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 20);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/Kova_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 20);
  vtuPath = "/lore/joshia5/Meshes/curved/Kova_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  fprintf(stderr, "start coarsening\n");
  for (LO adapt_itr = 0; adapt_itr < 7; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    coarsen_by_size(&mesh, opts);
  }
  //while ((coarsen_by_size(&mesh, opts)) && (mesh.nelems() > 50));
  fprintf(stderr, "finish coarsening\n");
  mesh.ask_qualities();
  writer = vtk::FullWriter("/lore/joshia5/Meshes/curved/kovaCoarsen_aft.vtk", &mesh);
  writer.write();
  fprintf(stderr, "finished kova case\n");
  
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
/*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/box_circleCut_4k_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/box_circleCut_4k_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/boxCircle_bef.vtk", &mesh);
  writer.write();
*/
 
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(100)));
  for (LO adapt_itr = 0; adapt_itr < 2; ++adapt_itr) 
    coarsen_by_size(&mesh, opts);
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

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_disc_collapse(&lib);
  //test_disc_validity(&lib);
  //test_tri_validity(&lib);
  //test_boxCircle_validity(&lib);
  //test_linear_tet_validity(&lib);
  //test_quadratic_tet_validity(&lib);
  //test_Kova_validity(&lib);
  //test_cubic_tet_validity(&lib);
  //test_collapse_kova(&lib);
  test_collapse_boxCircle(&lib);
  //test_collapse_cubicSlab(&lib);
  //test_antenna(&lib);

  return 0;
}
