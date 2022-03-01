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

void test_tet_validity(Library *lib) {
  auto mesh = Mesh(lib);
  auto comm = lib->world();
  binary::read(
      "/lore/joshia5/develop/mfem_omega/omega_h/meshes/Example_tet.osh",
      lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  auto dim = mesh.dim();
  auto coords = mesh.coords();
  auto edge_nCtrlPts = mesh.n_internal_ctrlPts(1);
  auto ev2v = mesh.ask_down(1, 0).ab2b;
  auto fe2e = mesh.ask_down(2, 1).ab2b;
  auto rv2v = mesh.ask_down(3, 0).ab2b;//0,2,1,3
  auto re2e = mesh.ask_down(3, 1).ab2b;//1,3,0,2,5,4
  auto rf2f = mesh.ask_down(3, 2).ab2b;//0,1,2,3

  auto vertCtrlPts = HostRead<Real>(Reals({0.0,0.0,0.0,
                                           1.0,0.0,0.0,
                                           0.0,1.0,0.0,
                                           0.0,0.0,1.0}));
  auto edgeCtrlPts = HostRead<Real>(Reals({
                                     2.0/3.0, 0.0, 0.0,
                                     1.0/3.0, 0.0, 0.0,
                                     //1
                                     1.0/3.0, 2.0/3.0, 0.0, 
                                     2.0/3.0, 1.0/3.0, 0.0,
                                     //3
                                     0.0, 1.0/3.0, 0.0, 
                                     0.0, 2.0/3.0, 0.0,
                                     //0
                                     0.0, 0.0, 1.0/3.0, 
                                     0.0, 0.0, 2.0/3.0,
                                     //2
                                     2.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 0.0, 2.0/3.0,
                                     //5
                                     0.0, 2.0/3.0, 1.0/3.0,
                                     0.0, 1.0/3.0, 2.0/3.0
                                     //4
                                     }));
  auto faceCtrlPts = HostRead<Real>(Reals({
                                     1.0/3.0, 1.0/3.0, 0.0,
                                     1.0/3.0, 0.0, 1.0/3.0,
                                     1.0/3.0, 1.0/3.0, 1.0/3.0,
                                     0.0, 1.0/3.0, 1.0/3.0}));

  Few<Real, 60> tet_pts;//ntet_pts*dim=20*3
  LO tet = 0;
  for (LO j = 0; j < vertCtrlPts.size(); ++j) {
    tet_pts[j] = vertCtrlPts[j];
  }
  for (LO j = 0; j < edgeCtrlPts.size(); ++j) {
    tet_pts[12 + j] = edgeCtrlPts[j];
  }
  for (LO j = 0; j < faceCtrlPts.size(); ++j) {
    tet_pts[48 + j] = faceCtrlPts[j];
  }
  Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);
  auto is_invalid = checkMinJacDet_3d(nodes_det);
  printf("tet is invalid %d %f %f %f %f\n", is_invalid, nodes_det[0],nodes_det[1],nodes_det[2],nodes_det[3]);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_disc_collapse(&lib);
  //test_disc_validity(&lib);
  //test_tri_validity(&lib);
  test_tet_validity(&lib);
  return 0;
}
