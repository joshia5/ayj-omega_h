#include <iostream>
#include <fstream>
#include <math.h>

#include<Omega_h_mesh.hpp>
#include<Omega_h_file.hpp>
#include<Omega_h_beziers.hpp>
#include<Omega_h_matrix.hpp>
#include<Omega_h_defines.hpp>
#include "Omega_h_for.hpp"

using namespace Omega_h;

Real B0(Real u) {
  return intpow(1-u, 3);
}

Real B1(Real u) {
  return 3*u*intpow(1-u, 2);
}

Real B2(Real u) {
  return 3*(1-u)*intpow(u, 2);
}

Real B3(Real u) {
  return intpow(u, 3);
}

void test_2d(Library *lib, const std::string &mesh_file, const char* vtu_file) {
  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_max_order(3);
  mesh.add_tags_for_ctrlPts();
  //make a pfor and populate the ctrl pts values in a write real sized
  //nents*numpts*dim
  //
  //set tags using the new arrays

  Real centerx = 0.5;
  Real centery = 0.5;
  Real R = 0.70711;

  // edge 1 (diameter, 0 to PI)
  Real xi_1 = 0.2748;
  Real xi_2 = 0.7252;
  Real s_c1 = PI * xi_1;
  Real s_c2 = PI * xi_2;

  Real x1 = centerx - R*cos(s_c1);
  Real y1 = centery + R*sin(s_c1);

  Real x2 = centerx - R*cos(s_c2);
  Real y2 = centery + R*sin(s_c2);

  Matrix<2,1> fx({x1, x2});
  Matrix<2,1> fy({y1, y2});

  Matrix<2,2> M1_inv ({B1(xi_1), B2(xi_1), B1(xi_2), B2(xi_2)});
  auto M1 = invert(M1_inv);
  Matrix<2,2> M2 ({B0(xi_1), B3(xi_1), B0(xi_2), B3(xi_2)});

  Matrix<2,1> p0({0, 1});//ask coords
  Matrix<2,1> p1({1, 0});
  auto Cx = M1*fx - M1*M2*p0;
  auto Cy = M1*fy - M1*M2*p1;

  auto u = Read<Real>(1000, 0.0, 0.001, "samplePts");
  auto u_h = HostRead<Real>(u);
  auto u_ex = Read<Real>(1000, 0.0 + PI/4, PI/1000, "exactPts");
  auto u_ex_h = HostRead<Real>(u_ex);
  std::ofstream edge1_file;
  edge1_file.open("edge1.csv");

  edge1_file << "x, y\n";
  for (LO i = 0; i < u.size(); ++i) {
    auto x_bezier = 0*B0(u[i]) + Cx(0,0)*B1(u[i]) + Cx(1,0)*B2(u[i]) + 1*B3(u[i]);
    auto y_bezier = 1*B0(u[i]) + Cy(0,0)*B1(u[i]) + Cy(1,0)*B2(u[i]) + 0*B3(u[i]);
    edge1_file << x_bezier << ", " << y_bezier << "\n";
  }
  edge1_file << "\n";
  edge1_file << "\n";
  edge1_file << "C1x, C1y\n";
  edge1_file << Cx(0,0) << ", " << Cy(0,0) << "\n";
  edge1_file << "C2x, C2y\n";
  edge1_file << Cx(1,0) << ", " << Cy(1,0) << "\n";
  edge1_file << "x1, y1\n";
  edge1_file << x1 << ", " << y1 << "\n";
  edge1_file << "x2, y2\n";
  edge1_file << x2 << ", " << y2 << "\n";
  edge1_file << "\n";
  edge1_file << "\n";
  edge1_file << "exact circle:\n";
  edge1_file << "x, y\n";
  for (LO i = 0; i < u_ex.size(); ++i) {
    auto x_circle = centerx - R*cos(u_ex[i]);
    auto y_circle = centery + R*sin(u_ex[i]);
    edge1_file << x_circle << ", " << y_circle << "\n";
  }

  edge1_file.close();
  //
  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  if (argc != 3) {
    Omega_h_fail("a.out <2d_in_osh> <2d_out_vtu>\n");
  };
  char const* path_2d = nullptr;
  char const* path_2d_vtu = nullptr;
  path_2d = argv[1];
  path_2d_vtu = argv[2];

  test_2d(&lib, path_2d, path_2d_vtu);

  return 0;
}
