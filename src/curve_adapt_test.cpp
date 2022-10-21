#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_curve_validity_3d.hpp>

using namespace Omega_h;

void test_adapt_inclusion(Library *lib) {
  auto comm = lib->world();

  auto mesh = binary::read(
      "/lore/joshia5/Meshes/curved/inclusion_3p_sizes.osh", comm);
                            
  for (LO i = 0; i <= mesh.dim(); ++i) {
    if (!mesh.has_tag(i, "global")) {
      mesh.add_tag(i, "global", 1, Omega_h::GOs(mesh.nents(i), 0, 1));
    }
  }

  vtk::FullWriter writer;

  /*
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
  */

  auto opts = AdaptOpts(&mesh);
  opts.should_swap = false;
  opts.should_refine = true;
  opts.should_filter_invalids = true;
  opts.verbosity = EXTRA_STATS;
  fprintf(stderr, "initial mesh size %d\n", mesh.nregions());
  I8 max_adapt_itr = 1;
  for (LO adapt_itr = 0; adapt_itr < max_adapt_itr; ++adapt_itr) {
    adapt(&mesh, opts);
  }
  auto qual = calc_crvQuality_3d(&mesh);
  /*
  writer = vtk::FullWriter(
      "../omega_h/meshes/boxCircle_aft.vtk", &mesh);
  writer.write();
  */
  return;
}

void test_cubic_tet_quality(Library *lib) {
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
  add_implied_metric_tag(&mesh);
  auto qual = HostRead<Real>(calc_crvQuality_3d(&mesh));
  printf("quality %f\n", qual[0]);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_adapt_inclusion(&lib);
  test_cubic_tet_quality(&lib);

  return 0;
}
