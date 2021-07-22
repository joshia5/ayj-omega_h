#include <iostream>

#include<Omega_h_mesh.hpp>
#include<Omega_h_file.hpp>
#include<Omega_h_beziers.hpp>

using namespace Omega_h;

void test_2d(Library *lib, const std::string &mesh_file, const char* vtu_file) {
  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);
  mesh.set_curved(1);
  mesh.set_order(3);

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
