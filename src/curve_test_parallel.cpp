#include <iostream>
#include <fstream>

#include<Omega_h_mesh.hpp>
#include<Omega_h_for.hpp>
#include<Omega_h_file.hpp>
#include<Omega_h_beziers.hpp>
#include<Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_array_ops.hpp>

using namespace Omega_h;

void test_sim_kova_quadratic(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/lore/joshia5/Meshes/curved/box_circleCut-100k_2p.osh", comm);
  vtk::FullWriter writer;
  writer = vtk::FullWriter("/lore/joshia5/Meshes/curved/box_circleCut-100k_2p_full.vtk", &mesh);
  writer.write();
  if (!mesh.has_tag(0, "bezier_pts"))
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);

  /*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe_3d(&mesh, &wireframe_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  */
  /*
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk_3d(&mesh, &curveVtk_mesh);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
*/
/*
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/KovaGeomSim-123tet_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh);
  vtuPath = "/lore/joshia5/Meshes/curved/KovaGeomSim-123tet_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
*/
  AdaptOpts opts(&mesh);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  auto desired_group_nelems = 1000000;
  while (nelems < desired_group_nelems) {
    if (!mesh.has_tag(0, "metric")) {
      add_implied_metric_tag(&mesh);
      adapt_refine(&mesh, opts);
      nelems = mesh.nglobal_ents(mesh.dim());
      std::cout << "mesh now has " << nelems << " total elements\n";
    }
    auto metrics = mesh.get_array<double>(0, "metric");
    metrics = multiply_each_by(metrics, 1.2);
    auto const metric_ncomps =
      divide_no_remainder(metrics.size(), mesh.nverts());
    mesh.add_tag(0, "metric", metric_ncomps, metrics);
    adapt_refine(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  /*
  vtk::write_parallel("/lore/joshia5/Meshes/curved/kova_refined.vtk", &mesh, 2);
  wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 2);
  vtuPath = "/lore/joshia5/Meshes/curved/KovaGeomSim-123tet_wireframe_refined.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  */
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 2);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/box_circleCut-100k_curvevtk_2p_"
    + std::to_string(comm->rank()) + ".vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  return;
}
/*
void test_disc(Library *lib) {
  auto comm = lib->world();
  //auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_2tri_order2.sms",
    //                        "/users/joshia5/Models/curved/disk_semi_geomsim.smd", comm);
  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
                            "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", comm);
  vtk::write_parallel("/lore/joshia5/Meshes/curved/disc_refine_100.vtk", &mesh, 2);

  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath = "/users/joshia5/Meshes/curved/disc100_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(comm->library());
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/users/joshia5/Meshes/curved/disc100_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  if (!mesh.has_tag(1, "global")) {
    mesh.add_tag(1, "global", 1, Omega_h::GOs(mesh.nedges(), 0, 1));
  }
  if (!mesh.has_tag(0, "global")) {
    mesh.add_tag(0, "global", 1, Omega_h::GOs(mesh.nverts(), 0, 1));
  }
  if (!mesh.has_tag(2, "global")) {
    mesh.add_tag(2, "global", 1, Omega_h::GOs(mesh.nfaces(), 0, 1));
  }
  AdaptOpts opts(&mesh);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  auto desired_group_nelems = 1000;
  while (nelems < desired_group_nelems) {
    if (!mesh.has_tag(0, "metric")) {
      add_implied_metric_tag(&mesh);
      adapt_refine(&mesh, opts);
      nelems = mesh.nglobal_ents(mesh.dim());
      std::cout << "mesh now has " << nelems << " total elements\n";
    }
    auto metrics = mesh.get_array<double>(0, "metric");
    metrics = multiply_each_by(metrics, 1.2);
    auto const metric_ncomps =
      divide_no_remainder(metrics.size(), mesh.nverts());
    mesh.add_tag(0, "metric", metric_ncomps, metrics);
    adapt_refine(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  vtk::write_parallel("/lore/joshia5/Meshes/curved/disc_refined_100to1k.vtk", &mesh, 2);
  wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/disc100to1k_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(comm->library());
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/disc100to1k_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  return;
}

*/
int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  /*
  if (argc != 6) {
    Omega_h_fail(
      "a.out <2d_in_osh> <2d_out_vtk> <3d_in_model-geomsim> <3d_in_mesh> <3d_out_vtk>\n");
  };*/

  //test_disc(&lib);
  test_sim_kova_quadratic(&lib);

  return 0;
}
