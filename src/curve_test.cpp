#include <iostream>

#include<Omega_h_mesh.hpp>
#include<Omega_h_file.hpp>
#include<Omega_h_beziers.hpp>
#include<Omega_h_matrix.hpp>
#include<Omega_h_defines.hpp>

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

  Matrix<2,1> p0({0, 1});
  Matrix<2,1> p1({1, 0});
  auto Cx = M1*fx - M1*M2*p0;
  auto Cy = M1*fy - M1*M2*p1;

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
