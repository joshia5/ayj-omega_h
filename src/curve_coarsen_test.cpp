#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_for.hpp>

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

  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_2tri_order2.sms",
                          "/users/joshia5/Models/curved/disk_semi_geomsim.smd", comm);
  //auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
    //                        "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", lib.world());
                            
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

  auto coords = mesh.coords();
  auto ev2v = mesh.get_adj(1, 0).ab2b;
  auto fe2e = mesh.get_adj(2, 1).ab2b;
  printf("coords %d, ev2v %d, fe2e %d\n", coords.size(), ev2v.size(), fe2e.size());
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(0.3)));
  while (coarsen_by_size(&mesh, opts));
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(0.6)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(),
        metric_eigenvalue_from_length(1.0)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.ask_qualities();
  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_disc_collapse(&lib);
  //test_disc_validity(&lib);
  //test_tri_validity(&lib);
  return 0;
}
