#include <Omega_h_file.hpp>
#include <Omega_h_filesystem.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_file(py::module& module) {
  py::class_<filesystem::path>(module, "path")
      .def(py::init<char const*>());
  //py::class_<Mesh>(module, "Mesh")
    //  .def(py::init<Mesh>());
  //filesystem::path (*path_c)(char const* source) = &filesystem::path::path;
  //filesystem::path (*path_s)(std::string const& source) = &filesystem::path::path; // typename not allowed
  Mesh (*gmsh_read_file)(filesystem::path const&, CommPtr) = &gmsh::read;
  void (*gmsh_write_file)(filesystem::path const&, Mesh*) = &gmsh::write;
  void (*vtk_write_vtu_dim)(std::string const&, Mesh*, Int, bool) =
      &vtk::write_vtu;
  void (*vtk_write_vtu)(std::string const&, Mesh*, bool) = &vtk::write_vtu;
  void (*vtk_write_parallel_dim)(std::string const&, Mesh*, Int, bool) =
      &vtk::write_parallel;
  void (*vtk_write_parallel)(std::string const&, Mesh*, bool) =
      &vtk::write_parallel;
  //Mesh (*binary_read_file)
    //(filesystem::path const&, CommPtr, bool) = &binary::read;
  Mesh (*binary_read_file)
    (std::string const&, CommPtr, bool) = &binary::read_s;
  void (*binary_write_file)
    (std::string const&, Mesh*) = &binary::write_s;
  module.def("gmsh_read_file", gmsh_read_file, "Read a Gmsh file");
  module.def("gmsh_write_file", gmsh_write_file, "Write a Gmsh file");
  module.def("vtk_write_vtu", vtk_write_vtu, "Write a mesh as a .vtu file",
      py::arg("path"), py::arg("mesh"), py::arg("compress") = true);
  module.def("vtk_write_vtu_dim", vtk_write_vtu_dim,
      "Write entities of one dimension as a .vtu file", py::arg("path"),
      py::arg("mesh"), py::arg("cell_dim"), py::arg("compress") = true);
  module.def("vtk_write_parallel", vtk_write_parallel,
      "Write a mesh as a directory of parallel VTK files", py::arg("path"),
      py::arg("mesh"), py::arg("compress") = true);
  module.def("vtk_write_parallel_dim", vtk_write_parallel_dim,
      "Write entities of one dimension as a directory of parallel VTK files",
      py::arg("path"), py::arg("mesh"), py::arg("cell_dim"),
      py::arg("compress") = true);
  module.def("binary_read_file", binary_read_file, "Read a .osh file",
      py::arg("path"), py::arg("comm"), py::arg("strict") = false);
  module.def("binary_write_file", binary_write_file, "Write a .osh file",
      py::arg("path"), py::arg("mesh"));
  //module.def("path_c", path_c, "Byte to path", py::arg("source"));
  //module.def("path_s", path_s, "String to path", py::arg("source"));
}

}  // namespace Omega_h
