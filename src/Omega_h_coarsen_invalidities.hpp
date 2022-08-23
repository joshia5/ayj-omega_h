#ifndef OMEGA_H_COARSEN_INVALIDITIES_HPP
#define OMEGA_H_COARSEN_INVALIDITIES_HPP

#include <Omega_h_mesh.hpp>

namespace Omega_h {

class Mesh;

LOs coarsen_invalidities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

LOs coarsen_invalidities_new_mesh(Mesh* mesh, LOs cands2edges,
  Read<I8> cand_codes, Mesh* new_mesh, LOs const old_verts2new_verts,
  Read<I8> const verts_are_keys, LOs const keys2verts_onto,
  LOs const prods2new_ents);

}  // end namespace Omega_h

#endif
