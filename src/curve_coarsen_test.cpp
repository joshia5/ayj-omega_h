#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>

#include <cmath>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  auto mesh = meshsim::read("/users/joshia5/Meshes/curved/disk_semi_100_order2.sms",
                            "/users/joshia5/Models/curved/disk_semi_geomsim_100.smd", lib.world());
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);

  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.3)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.6)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.0)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.ask_qualities();
  Omega_h::vtk::write_parallel("disc_8tri_coarsen.vtk", &mesh, mesh.dim());
  return 0;
}
