#include "Omega_h_collapse.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

void make_cavity_class(Mesh* mesh, LOs const keys2verts) {
  auto nkeys = keys2verts.size();
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto faces2class_dim = mesh->get_array<I8>(FACE, "class_dim");
  auto v2vt = (mesh->ask_up(0, 3)).a2ab;//t tet
  auto vt2t = (mesh->ask_up(0, 3)).ab2b;
  auto t2tf = (mesh->get_adj(3, 2)).ab2b;
  auto f2fe = (mesh->get_adj(2, 1)).ab2b;
  auto e2ev = (mesh->get_adj(1, 0)).ab2b;
  auto f2ft = (mesh->ask_up(2, 3)).a2ab;
  auto ft2t = (mesh->ask_up(2, 3)).ab2b;

  auto cav_classdim_f = Write<I8>(mesh->nfaces());
  auto cav_classdim_e = Write<I8>(mesh->nedges());
  auto cav_classdim_v = Write<I8>(mesh->nverts());

  auto cav_tets = Write<I8>(mesh->nregions());
  auto f = OMEGA_H_LAMBDA(LO key) {
    LO const v_key= keys2verts[key];
    for (LO vt = v2vt[v_key]; vt < v2vt[v_key+1]; ++vt) {
      //adj tets of vkey
      LO const adj_t = vt2t[vt];
      //mark all cavity tets
      cav_tets[adj_t] = 1;
    }
    for (LO vt = v2vt[v_key]; vt < v2vt[v_key+1]; ++vt) {
      LO const adj_t = vt2t[vt];
      //itr over faces
      for (LO f=0; f<4; ++f) {
        LO const down_face = t2tf[adj_t*4+f];
        LO n_adj_cav_tets = 0;
        for (LO ft=f2ft[down_face]; ft<f2ft[down_face+1]; ++ft) {
          LO const adj_t2 = ft2t[ft];
          if (cav_tets[adj_t2] == 1) {
            ++n_adj_cav_tets;
          }
        }
        if (n_adj_cav_tets == 1) {
          fprintf(stderr, "bdry face %d\n", down_face);
          cav_classdim_f[down_face] = 2;
          // if adj tets in cavity = 1;
          // face is on bdry
          for (LO e=0; e<3; ++e) {
            LO const down_edge = f2fe[down_face*3+e];
            cav_classdim_e[down_edge] = 1;
            // mark down adj edge and verts of this face as being on bdry
            for (LO v=0; v<2; ++v) {
              LO const down_vert = e2ev[down_edge*2+v];
              cav_classdim_v[down_vert] = 1;
            }
          }
        }
      }
    }
  };
  parallel_for(nkeys, f, "make_cavity_class");
  mesh->add_tag<I8>(0, "cav_classdim", 1, Read<I8>(cav_classdim_v));
  mesh->add_tag<I8>(1, "cav_classdim", 1, Read<I8>(cav_classdim_e));
  mesh->add_tag<I8>(2, "cav_classdim", 1, Read<I8>(cav_classdim_f));

  return;
}

}  // end namespace Omega_h
