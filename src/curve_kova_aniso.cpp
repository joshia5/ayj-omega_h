#include <iostream>
#include <fstream>
#include <math.h>

#include <Omega_h_timer.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_refine.hpp>
#include <Omega_h_array_ops.hpp>

using namespace Omega_h;

void test_sim_kova_quadratic(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet.osh", comm);
  if (!mesh.has_tag(0, "bezier_pts")) 
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);

  elevate_curve_order_2to3(&mesh);
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }
  AdaptOpts opts(&mesh);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  auto desired_group_nelems = 5000;
  Now t0 = now();
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
    refine_by_size(&mesh, opts);
    //adapt_refine(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  Now t1 = now();
  std::cout << "total refine time: " << (t1 - t0) << " seconds\n";
  //vtk::write_parallel("/lore/joshia5/Meshes/curved/kova_refined.vtk", &mesh, 2);
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 10);
  std::string vtuPath = 
    "/lore/joshia5/Meshes/curved/Kova-123tetwire_ansiorefin.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 10);
  vtuPath = 
    "/lore/joshia5/Meshes/curved/Kova-123tetcurveVtk_anisorefin.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);
  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  /*
  if (argc != 6) {
    Omega_h_fail(
      "a.out <2d_in_osh> <2d_out_vtk> <3d_in_model-geomsim> <3d_in_mesh> <3d_out_vtk>\n");
  };
  char const* path_2d = nullptr;
  char const* path_2d_vtk = nullptr;
  path_2d = argv[1];
  path_2d_vtk = argv[2];

  char const* path_3d_g = nullptr;
  char const* path_3d_m = nullptr;
  char const* path_3d_vtk = nullptr;
  path_3d_g = argv[3];
  path_3d_m = argv[4];
  path_3d_vtk = argv[5];

  */
  test_sim_kova_quadratic(&lib);

  return 0;
}
