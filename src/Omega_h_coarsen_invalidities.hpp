#ifndef OMEGA_H_COARSEN_INVALIDITIES_HPP
#define OMEGA_H_COARSEN_INVALIDITIES_HPP

#include <Omega_h_mesh.hpp>

namespace Omega_h {

class Mesh;

LOs coarsen_invalidities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

LOs coarsen_invalidities_new_ops(Mesh* mesh, LOs cands2edges,
  Read<I8> cand_codes, Mesh* new_mesh, LOs const old_verts2new_verts);

}  // end namespace Omega_h

#endif
