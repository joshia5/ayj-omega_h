#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"

namespace Omega_h {

template <Int dim>
void coarsen_curved_verts_and_edges(Mesh *mesh, Mesh *new_mesh, LOs old2new,
    LOs prods2new, LOs old_verts2new_verts, LOs keys2verts, LOs keys2verts_onto,
    LOs keys2prods) {
  auto const nold_verts = mesh->nverts();
  auto const nold_edges = mesh->nedges();
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  if (!mesh->has_tag(0, "bezier_pts")) {
    mesh->add_tag<Real>(0, "bezier_pts", dim, mesh->coords());
  }
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<Real> vert_ctrlPts(nnew_verts*1*dim, INT8_MAX);

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
  new_mesh->add_tag<Real>(0, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(0, Reals(vert_ctrlPts));

  //copy ctrl pts for edges
  auto copy_sameedgePts = OMEGA_H_LAMBDA(LO i) {
    if (old2new[i] != -1) {
      LO new_edge = old2new[i];
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*n_edge_pts*dim + d] = old_edgeCtrlPts[i*n_edge_pts*dim + d];
        edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = old_edgeCtrlPts[i*n_edge_pts*dim + dim + d];
      }
    }
  };
  parallel_for(nold_edges, std::move(copy_sameedgePts),
      "copy_same_edgectrlPts");
  auto prod_edge_points = OMEGA_H_LAMBDA(LO i) {
    LO e = prods2new[i];
    auto v0 = new_ev2v[e*2 + 0];
    auto v1 = new_ev2v[e*2 + 1];
    for (LO j=0; j<dim; ++j) {
      edge_ctrlPts[e*n_edge_pts*dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*(1.0/3.0);
          //TODO verify that new pts should be equidistant
          //(new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_1_cube();
      edge_ctrlPts[e*n_edge_pts*dim + dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*(2.0/3.0);
          //(new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_2_cube();
    }
  };
  parallel_for(prods2new.size(), std::move(prod_edge_points),
      "prod_edge_points");

  auto newvert_gdim =  new_mesh->get_array<I8>(0, "class_dim");
  auto newedge_gdim =  new_mesh->get_array<I8>(1, "class_dim");
  auto newedge_gid =  new_mesh->get_array<LO>(1, "class_id");
  auto oldvert_gdim =  mesh->get_array<I8>(0, "class_dim");
  auto oldvert_gid =  mesh->get_array<LO>(0, "class_id");
  auto oldedge_gid =  mesh->get_array<LO>(1, "class_id");

  Write<LO> new_edge2keys(new_mesh->nedges(), -1);
  Write<LO> new_edge2mid_pt(new_mesh->nedges()*dim, -1.0);

  auto v2v_old = mesh->ask_star(0);
  auto v2vv_old = v2v_old.a2ab;
  auto vv2v_old = v2v_old.ab2b;
  auto v2v_new = new_mesh->ask_star(0);
  auto v2vv_new = v2v_new.a2ab;
  auto vv2v_new = v2v_new.ab2b;
  /*
  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   new_mesh->nverts());
  auto nv2ov_ab2b = new_verts2old_verts.ab2b;
  auto nv2ov_a2ab = new_verts2old_verts.a2ab;
  */
  auto find_bdry_edges = OMEGA_H_LAMBDA(LO i) {
    LO v_key = keys2verts[i];
    LO v_onto = keys2verts_onto[i];
    if ((oldvert_gdim[v_key] == 1) && (oldvert_gdim[v_onto] == 1) && 
        (oldvert_gid[v_key] == oldvert_gid[v_onto])) {
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO new_edge = prods2new[prod];
        if ((newedge_gdim[new_edge] == 1) && (newedge_gid[new_edge] == oldvert_gid[v_key])) {
          new_edge2keys[new_edge] = v_key;
          break;
        } 
      }

    }
      /*
    for (LO vv = v2vv_old[v_key]; vv < v2vv_old[v_key+1]; ++vv) {
      LO adj_v = vv2v_old[vv];
      for (LO vv2 = v2vv_new[v_onto_new]; vv2 < v2vv_new[v_onto_new+1]; ++vv2) {
        LO adj_v2_new = vv2v_new[vv2];
        OMEGA_H_CHECK((nv2ov_a2ab[adj_v2_new+1] - nv2ov_a2ab[adj_v2_new]) == 1);
        LO adj_v2 = nv2ov_ab2b[nv2ov_a2ab[adj_v2_new]];
        if (adj_v2 == adj_v) break;
      }
    }
    */
  };
  parallel_for(keys2verts.size(), std::move(find_bdry_edges));

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));

  return;
}

template <Int dim>
void coarsen_curved_faces(Mesh *mesh, Mesh *new_mesh, LOs old2new,
    LOs prods2new) {
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_faces = new_mesh->nfaces();
  auto const nold_faces = mesh->nfaces();
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  Write<Real> face_ctrlPts(nnew_faces*1*dim, INT8_MAX);

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
      "copy_same_facectrlPts");

  auto face_centroids = OMEGA_H_LAMBDA(LO i) {
    LO f = prods2new[i];
    auto v0 = new_fv2v[f*3 + 0];
    auto v1 = new_fv2v[f*3 + 1];
    auto v2 = new_fv2v[f*3 + 2];
    for (LO j=0; j<dim; ++j) {
      face_ctrlPts[f*dim + j] = (new_coords[v0*dim + j] +
          new_coords[v1*dim + j] + new_coords[v2*dim + j])/3.0;
    }
  };
  parallel_for(prods2new.size(), std::move(face_centroids),
      "face_centroids");
  new_mesh->add_tag<Real>(2, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(2, Reals(face_ctrlPts));

  return;
}

#define OMEGA_H_INST(T)                                                       
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h
