#include <iostream>
#include <fstream>

#include <Omega_h_timer.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_array_ops.hpp>

using namespace Omega_h;

void test_sim_kova_quadratic(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/lore/joshia5/develop/mfem_omega/omega_h/meshes/box_circleCut-3p5mil_2p.osh", comm);
  if (!mesh.has_tag(0, "bezier_pts"))
    mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);

  AdaptOpts opts(&mesh);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  auto desired_group_nelems = 10000000;
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
    adapt_refine(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  Now t1 = now();
  std::cout << "total refine time: " << (t1 - t0) << " seconds\n";
  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &curveVtk_mesh, 2);
  std::string vtuPath = "/lore/joshia5/Meshes/curved/box_circleCut-100k_curvevtk_2p_"
    + std::to_string(comm->rank()) + ".vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_sim_kova_quadratic(&lib);

  return 0;
}
