#include "Omega_h_collapse.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

void make_cavity_class(Mesh* mesh, LOs keys2verts) {
  auto nkeys = keys2verts.size();
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto faces2class_dim = mesh->get_array<I8>(FACE, "class_dim");
  auto v2vt = (mesh->ask_up(0, 3)).a2ab;
  auto v2vt = (mesh->ask_up(0, 3)).ab2b;

  auto cav_classdim_f = Write<I8>(mesh->nfaces());
  auto cav_classdim_e = Write<I8>(mesh->nedges());
  auto cav_classdim_v = Write<I8>(mesh->nverts());

  auto f = OMEGA_H_LAMBDA(LO key) {
    LO const vert = keys2verts[key];
    for (LO vt = v2vt[v_key]; vt < v2vt[v_key + 1]; ++vt) {
      //adj tets of vkey
      LO adj_t = vt2t[vt];
    }
    //mark all cavity tets
    //itr over faces
      // if adj tets in cavity = 1; 
        // face is on bdry
        // mark down adj edge and verts of this face as being on bdry
  };
  parallel_for(nkeys, f, "make_cavity_class");
}

}  // end namespace Omega_h
