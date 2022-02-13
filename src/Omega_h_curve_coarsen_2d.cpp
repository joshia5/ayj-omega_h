#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_matrix.hpp"

namespace Omega_h {

LOs coarsen_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs old_verts2new_verts,
                                     LOs old_edges2new_edges) {
  printf("in coarsen curved edges fn\n");
  OMEGA_H_TIME_FUNCTION;
  //auto const nold_edge = old2new.size();
  auto const nold_verts = mesh->nverts();
  auto const nold_faces = mesh->nfaces();
  auto const nold_edges = mesh->nedges();
  OMEGA_H_CHECK(nold_verts == old_verts2new_verts.size());
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  if (!mesh->has_tag(0, "bezier_pts")) {
    mesh->add_tag<Real>(0, "bezier_pts", mesh->dim(), mesh->coords());
  }
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  //auto const order = mesh->get_max_order();
  auto const dim = mesh->dim();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const nnew_verts = new_mesh->nverts();
  auto const nnew_face = new_mesh->nfaces();
  auto const n_face_pts = mesh->n_internal_ctrlPts(2);
  auto new_fv2v = mesh->ask_down(2, 0).ab2b;

  Write<Real> face_ctrlPts(nnew_face*n_face_pts*dim, INT8_MAX);
  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<Real> vert_ctrlPts(nnew_verts*1*dim, INT8_MAX);

  //auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
    //                                               nnew_verts);

  printf("curve coarsen L44, nkeys %d, old2new.size %d\n",
      keys2prods.size(), old2new.size());
  //copy ctrl pts for verts
  auto copy_sameCtrlPts = OMEGA_H_LAMBDA(LO i) {
    if (old_verts2new_verts[i] != -1) {
      LO new_vert = old_verts2new_verts[i];
      for (I8 d = 0; d < dim; ++d) {
        vert_ctrlPts[new_vert*dim + d] = old_vertCtrlPts[i*dim + d];
      }
    }
  };
  parallel_for(nold_verts, std::move(copy_sameCtrlPts),
      "copy same vtx ctrlPts");
  new_mesh->add_tag<Real>(0, "bezier_pts", mesh->dim());
  new_mesh->set_tag_for_ctrlPts(0, Reals(vert_ctrlPts));

  //copy ctrl pts for faces
  auto copy_samefacePts = OMEGA_H_LAMBDA(LO i) {
    if (old2new[i] != -1) {
      LO new_face = old2new[i];
      for (I8 d = 0; d < dim; ++d) {
        face_ctrlPts[new_face*dim + d] = old_faceCtrlPts[i*dim + d];
      }
    }
  };
  parallel_for(nold_faces, std::move(copy_samefacePts),
      "copy same facectrlPts");
  auto face_centroids = OMEGA_H_LAMBDA(LO i) {
    auto v0 = new_fv2v[i*3 + 0];
    auto v1 = new_fv2v[i*3 + 1];
    auto v2 = new_fv2v[i*3 + 2];
    for (LO j=0; j<dim; ++j) {
      face_ctrlPts[i*dim + j] = (new_coords[v0*dim + j] +
          new_coords[v1*dim + j] + new_coords[v2*dim + j])/3.0;
    }
  };
  parallel_for(nnew_face, std::move(face_centroids),
      "face_centroids");
  new_mesh->add_tag<Real>(2, "bezier_pts", n_face_pts*dim);
  new_mesh->set_tag_for_ctrlPts(2, Reals(face_ctrlPts));

  auto copy_sameedgePts = OMEGA_H_LAMBDA(LO i) {
    if (old_edges2new_edges[i] != -1) {
      LO new_edge = old_edges2new_edges[i];
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*dim + d] = old_edgeCtrlPts[i*dim + d];
      }
    }
  };
  parallel_for(nold_edges, std::move(copy_sameedgePts),
      "copy same edgectrlPts");

  //TODO change this to work on only prod edges
  auto prod_edge_points = OMEGA_H_LAMBDA(LO i) {
    auto v0 = new_ev2v[i*2 + 0];
    auto v1 = new_ev2v[i*2 + 1];
    for (LO j=0; j<dim; ++j) {
      edge_ctrlPts[i*n_edge_pts*dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_1_cube();
      edge_ctrlPts[i*n_edge_pts*dim + dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_2_cube();
    }
  };
  parallel_for(nnew_edge, std::move(prod_edge_points),
      "edge_points");
  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));

  auto valid_tris = checkValidity(new_mesh, prods2new, 2);
  fprintf(stderr, "out of check validity\n");

  return LOs();
}

#define OMEGA_H_INST(T)
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h
