#include <iostream>
#include <fstream>
#include <math.h>

#include<Omega_h_mesh.hpp>
#include<Omega_h_file.hpp>
#include<Omega_h_beziers.hpp>
#include<Omega_h_matrix.hpp>
#include<Omega_h_defines.hpp>
#include<Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_array_ops.hpp>

using namespace Omega_h;

void test_linearTri_toCubicCircle(Library *lib, const std::string &mesh_file,
                                  const char* vtk_file) {
  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  auto coords = mesh.coords();
  auto coords_h = HostRead<Real>(coords);
  auto edge_nCtrlPts = mesh.n_internal_ctrlPts(1);
  auto dim = mesh.dim();
  HostWrite<Real> edge_ctrlPts(mesh.nedges()*edge_nCtrlPts*dim);
  auto ev2v = mesh.ask_down(1, 0).ab2b;
  
  Real centerx = 0.5;
  Real centery = 0.5;
  Real R = 0.70711;
  Real xi_1 = 0.2748;
  Real xi_2 = 0.7252;

  std::ofstream edge_file;
  edge_file.open("edges.csv");

  // edge 1(1) (diameter, starts at pi/4, on for pi)
  {
    Real s_c1 = PI * xi_1;
    Real s_c2 = PI * xi_2;
    Matrix<2,1> p0({coords_h[4], coords_h[5]});
    Matrix<2,1> p1({coords_h[2], coords_h[3]});
    Real x1 = centerx - R*cos(PI/4 + s_c1);
    Real y1 = centery + R*sin(PI/4 + s_c1);
    Real x2 = centerx - R*cos(PI/4 + s_c2);
    Real y2 = centery + R*sin(PI/4 + s_c2);

    Matrix<2,1> fx({x1, x2});
    Matrix<2,1> fy({y1, y2});

    Matrix<2,2> M1_inv ({B1_cube(xi_1), B2_cube(xi_1), B1_cube(xi_2),
                         B2_cube(xi_2)});
    auto M1 = invert(M1_inv);
    Matrix<2,2> M2 ({B0_cube(xi_1), B3_cube(xi_1), B0_cube(xi_2),
                     B3_cube(xi_2)});

    auto Cx = M1*fx - M1*M2*p0;
    auto Cy = M1*fy - M1*M2*p1;

    edge_ctrlPts[1*edge_nCtrlPts*dim + 0*dim + 0] = Cx(0,0);
    edge_ctrlPts[1*edge_nCtrlPts*dim + 0*dim + 1] = Cy(0,0);
    edge_ctrlPts[1*edge_nCtrlPts*dim + 1*dim + 0] = Cx(1,0);
    edge_ctrlPts[1*edge_nCtrlPts*dim + 1*dim + 1] = Cy(1,0);

    auto u = Read<Real>(1000, 0.0, 0.001, "samplePts");
    auto u_h = HostRead<Real>(u);
    auto u_ex = Read<Real>(1000, PI/4, 2*PI/1000, "exactPts");
    auto u_ex_h = HostRead<Real>(u_ex);

    edge_file << "x, y\n";
    for (LO i = 0; i < u.size(); ++i) {
      auto x_bezier = p0(0,0)*B0_cube(u_h[i]) + Cx(0,0)*B1_cube(u_h[i]) + 
                      Cx(1,0)*B2_cube(u_h[i]) + p1(0,0)*B3_cube(u_h[i]);
      auto y_bezier = p0(1,0)*B0_cube(u_h[i]) + Cy(0,0)*B1_cube(u_h[i]) +
                      Cy(1,0)*B2_cube(u_h[i]) + p1(1,0)*B3_cube(u_h[i]);
      edge_file << x_bezier << ", " << y_bezier << "\n";
    }
    edge_file << "\n";
    edge_file << "\n";
    edge_file << "C1x, C1y\n";
    edge_file << Cx(0,0) << ", " << Cy(0,0) << "\n";
    edge_file << "C2x, C2y\n";
    edge_file << Cx(1,0) << ", " << Cy(1,0) << "\n";
    edge_file << "x1, y1\n";
    edge_file << x1 << ", " << y1 << "\n";
    edge_file << "x2, y2\n";
    edge_file << x2 << ", " << y2 << "\n";
    edge_file << "\n";
    edge_file << "\n";
    edge_file << "exact circle:\n";
    edge_file << "x, y\n";
    for (LO i = 0; i < u_ex.size(); ++i) {
      auto x_circle = centerx - R*cos(u_ex_h[i]);
      auto y_circle = centery + R*sin(u_ex_h[i]);
      edge_file << x_circle << ", " << y_circle << "\n";
    }
  }
  edge_file << "\n";

  //

  // edge 2(0) (along x, starts at 5*pi/4 goes on for pi/2)
  {
    Real s_c1 = PI/2 * xi_1;
    Real s_c2 = PI/2 * xi_2;
    Matrix<2,1> p0({coords_h[2], coords_h[3]});
    Matrix<2,1> p1({coords_h[0], coords_h[1]});
    Real x1 = centerx - R*cos(5*PI/4 + s_c1);
    Real y1 = centery + R*sin(5*PI/4 + s_c1);
    Real x2 = centerx - R*cos(5*PI/4 + s_c2);
    Real y2 = centery + R*sin(5*PI/4 + s_c2);

    Matrix<2,1> fx({x1, x2});
    Matrix<2,1> fy({y1, y2});

    Matrix<2,2> M1_inv ({B1_cube(xi_1), B2_cube(xi_1), B1_cube(xi_2), B2_cube(xi_2)});
    auto M1 = invert(M1_inv);
    Matrix<2,2> M2 ({B0_cube(xi_1), B3_cube(xi_1), B0_cube(xi_2), B3_cube(xi_2)});

    auto Cx = M1*fx - M1*M2*p0;
    auto Cy = M1*fy - M1*M2*p1;

    edge_ctrlPts[0*edge_nCtrlPts*dim + 0*dim + 0] = Cx(0,0);
    edge_ctrlPts[0*edge_nCtrlPts*dim + 0*dim + 1] = Cy(0,0);
    edge_ctrlPts[0*edge_nCtrlPts*dim + 1*dim + 0] = Cx(1,0);
    edge_ctrlPts[0*edge_nCtrlPts*dim + 1*dim + 1] = Cy(1,0);

    auto u = Read<Real>(1000, 0.0, 0.001, "samplePts");
    auto u_h = HostRead<Real>(u);

    edge_file << "for edge 2\n";
    edge_file << "x, y\n";
    for (LO i = 0; i < u.size(); ++i) {
      auto x_bezier = p0(0,0)*B0_cube(u_h[i]) + Cx(0,0)*B1_cube(u_h[i]) +
                      Cx(1,0)*B2_cube(u_h[i]) + p1(0,0)*B3_cube(u_h[i]);
      auto y_bezier = p0(1,0)*B0_cube(u_h[i]) + Cy(0,0)*B1_cube(u_h[i]) +
                      Cy(1,0)*B2_cube(u_h[i]) + p1(1,0)*B3_cube(u_h[i]);
      edge_file << x_bezier << ", " << y_bezier << "\n";
    }
    edge_file << "\n";
    edge_file << "\n";
    edge_file << "C1x, C1y\n";
    edge_file << Cx(0,0) << ", " << Cy(0,0) << "\n";
    edge_file << "C2x, C2y\n";
    edge_file << Cx(1,0) << ", " << Cy(1,0) << "\n";
    edge_file << "x1, y1\n";
    edge_file << x1 << ", " << y1 << "\n";
    edge_file << "x2, y2\n";
    edge_file << x2 << ", " << y2 << "\n";
  }
  edge_file << "\n";

  // edge 3(2) (along y, starts at 7*pi/4 goes on for pi/2)
  {
    Real s_c1 = PI/2 * xi_1;
    Real s_c2 = PI/2 * xi_2;
    Matrix<2,1> p0({coords_h[0], coords_h[1]});
    Matrix<2,1> p1({coords_h[4], coords_h[5]});
    Real x1 = centerx - R*cos(7*PI/4 + s_c1);
    Real y1 = centery + R*sin(7*PI/4 + s_c1);
    Real x2 = centerx - R*cos(7*PI/4 + s_c2);
    Real y2 = centery + R*sin(7*PI/4 + s_c2);

    Matrix<2,1> fx({x1, x2});
    Matrix<2,1> fy({y1, y2});

    Matrix<2,2> M1_inv ({B1_cube(xi_1), B2_cube(xi_1), B1_cube(xi_2),
                         B2_cube(xi_2)});
    auto M1 = invert(M1_inv);
    Matrix<2,2> M2 ({B0_cube(xi_1), B3_cube(xi_1), B0_cube(xi_2),
                     B3_cube(xi_2)});

    auto Cx = M1*fx - M1*M2*p0;
    auto Cy = M1*fy - M1*M2*p1;

    edge_ctrlPts[2*edge_nCtrlPts*dim + 0*dim + 0] = Cx(0,0);
    edge_ctrlPts[2*edge_nCtrlPts*dim + 0*dim + 1] = Cy(0,0);
    edge_ctrlPts[2*edge_nCtrlPts*dim + 1*dim + 0] = Cx(1,0);
    edge_ctrlPts[2*edge_nCtrlPts*dim + 1*dim + 1] = Cy(1,0);

    auto u = Read<Real>(1000, 0.0, 0.001, "samplePts");
    auto u_h = HostRead<Real>(u);

    edge_file << "for edge 3\n";
    edge_file << "x, y\n";
    for (LO i = 0; i < u.size(); ++i) {
      auto x_bezier = p0(0,0)*B0_cube(u_h[i]) + Cx(0,0)*B1_cube(u_h[i]) +
                      Cx(1,0)*B2_cube(u_h[i]) + p1(0,0)*B3_cube(u_h[i]);
      auto y_bezier = p0(1,0)*B0_cube(u_h[i]) + Cy(0,0)*B1_cube(u_h[i]) +
                      Cy(1,0)*B2_cube(u_h[i]) + p1(1,0)*B3_cube(u_h[i]);
      edge_file << x_bezier << ", " << y_bezier << "\n";
    }
    edge_file << "\n";
    edge_file << "\n";
    edge_file << "C1x, C1y\n";
    edge_file << Cx(0,0) << ", " << Cy(0,0) << "\n";
    edge_file << "C2x, C2y\n";
    edge_file << Cx(1,0) << ", " << Cy(1,0) << "\n";
    edge_file << "x1, y1\n";
    edge_file << x1 << ", " << y1 << "\n";
    edge_file << "x2, y2\n";
    edge_file << x2 << ", " << y2 << "\n";
  }

  edge_file.close();

  mesh.set_tag_for_ctrlPts(1, Reals(edge_ctrlPts.write()));

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
  auto desired_group_nelems = 2;
  while (nelems < desired_group_nelems) {
    if (!mesh.has_tag(0, "metric")) {
      add_implied_metric_tag(&mesh);
      adapt(&mesh, opts);
      nelems = mesh.nglobal_ents(mesh.dim());
      std::cout << "mesh now has " << nelems << " total elements\n";
    }
    auto metrics = mesh.get_array<double>(0, "metric");
    metrics = multiply_each_by(metrics, 1.2);
    auto const metric_ncomps =
      divide_no_remainder(metrics.size(), mesh.nverts());
    mesh.add_tag(0, "metric", metric_ncomps, metrics);
    adapt(&mesh, opts);
    nelems = mesh.nglobal_ents(mesh.dim());
    std::cout << "mesh now has " << nelems << " total elements\n";
  }
  vtk::write_parallel("/users/joshia5/Meshes/curved/4tri.vtk", &mesh, 2);

  return;
}

void test_sim_linearToCubic(Library *lib, const std::string &model_file,
                            const std::string &mesh_file,
                            const char* vtk_file) {
  auto comm = lib->world();
  auto mesh = meshsim::read(mesh_file, model_file, comm);
  calc_quad_ctrlPts_from_interpPts(&mesh);
 
  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe(&mesh, &wireframe_mesh);
  std::string vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk(&mesh, &curveVtk_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  elevate_curve_order_2to3(&mesh);

  auto cubic_wireframe_mesh = Mesh(comm->library());
  cubic_wireframe_mesh.set_comm(comm);
  build_cubic_wireframe(&mesh, &cubic_wireframe_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_cubic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_wireframe_mesh, 1);

  auto cubic_curveVtk_mesh = Mesh(comm->library());
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk(&mesh, &cubic_curveVtk_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_cubic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  elevate_curve_order_3to4(&mesh);

  auto quartic_wireframe_mesh = Mesh(comm->library());
  quartic_wireframe_mesh.set_comm(comm);
  build_quartic_wireframe(&mesh, &quartic_wireframe_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_quartic_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &quartic_wireframe_mesh, 1);

  auto quartic_curveVtk_mesh = Mesh(comm->library());
  quartic_curveVtk_mesh.set_comm(comm);
  build_quartic_curveVtk(&mesh, &quartic_curveVtk_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/box_circleCut-30reg_quartic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &quartic_curveVtk_mesh, 2);

  elevate_curve_order_4to5(&mesh);
  elevate_curve_order_5to6(&mesh);

  return;
}

void test_sim_kova_quadratic(Library *lib) {
  auto comm = lib->world();
  auto mesh = binary::read("/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet.osh", comm);
  calc_quad_ctrlPts_from_interpPts(&mesh);

  auto wireframe_mesh = Mesh(comm->library());
  wireframe_mesh.set_comm(comm);
  build_quadratic_wireframe(&mesh, &wireframe_mesh);
  std::string vtuPath = "/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet_wireframe.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);

  auto curveVtk_mesh = Mesh(comm->library());
  curveVtk_mesh.set_comm(comm);
  build_quadratic_curveVtk(&mesh, &curveVtk_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/KovaGeomSim-quadratic_123tet_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &curveVtk_mesh, 2);

  elevate_curve_order_2to3(&mesh);
  elevate_curve_order_3to4(&mesh);
  auto quartic_curveVtk_mesh = Mesh(lib);
  quartic_curveVtk_mesh.set_comm(comm);
  build_quartic_curveVtk(&mesh, &quartic_curveVtk_mesh);
  vtuPath = "/users/joshia5/Meshes/curved/KovaGeomSim_quartic_curveVtk.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &quartic_curveVtk_mesh, 2);

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

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

  test_linearTri_toCubicCircle(&lib, path_2d, path_2d_vtk);
  test_sim_linearToCubic(&lib, path_3d_g, path_3d_m, path_3d_vtk);
  test_sim_kova_quadratic(&lib);

  return 0;
}
