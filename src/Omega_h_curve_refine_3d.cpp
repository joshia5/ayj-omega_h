#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_matrix.hpp"
#include "Omega_h_atomics.hpp"

namespace Omega_h {

  /*
LO max_degree_key2oldents(LO const nold_edge, LOs old_key2keyent, LOs keys2edges) {
  auto keyedges_noldents_w = Write<LO>(nold_edge, -1);
  auto count_oldents = OMEGA_H_LAMBDA (LO key) {
    if (old2new[old_edge] == -1) {//old edge is key
      //get num up adj ents
      LO const num_adj_ents = old_key2keyent[old_edge + 1] - old_key2keyent[old_edge];
      keyedges_noldents_w[old_edge] = num_adj_ents;
    }
  };
  parallel_for(keys2edges.size(), std::move(count_oldents));
  return get_max(LOs(keyedges_noldents_w));
}
*/

LOs create_curved_verts_and_edges_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges) {
  auto const order = mesh->get_max_order();
  OMEGA_H_CHECK(order >= 2);
  OMEGA_H_TIME_FUNCTION;
  auto const nold_edge = old2new.size();
  auto const nold_verts = mesh->nverts();
  OMEGA_H_CHECK(nold_verts == old_verts2new_verts.size());
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const dim = mesh->dim();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim);
  Write<Real> vert_ctrlPts(nnew_verts*1*dim);

  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  LO max_degree_key2oldface = -1;
  {
    auto keyedges_noldfaces_w = Write<LO>(nold_edge, -1);
    auto count_oldfaces = OMEGA_H_LAMBDA (LO old_edge) {
      if (old2new[old_edge] == -1) {
        LO const num_adj_faces = old_e2ef[old_edge + 1] - old_e2ef[old_edge];
        keyedges_noldfaces_w[old_edge] = num_adj_faces;
      }
    };
    parallel_for(nold_edge, std::move(count_oldfaces));
    max_degree_key2oldface = get_max(LOs(keyedges_noldfaces_w));
  }

  auto nkeys = keys2edges.size();
  OMEGA_H_CHECK(order == 3);
  auto keys2old_faces_w = Write<LO>(max_degree_key2oldface*nkeys, -1);

  auto create_crv_prod_edges = OMEGA_H_LAMBDA (LO key) {
    LO const old_edge = keys2edges[key];
    LO const v0_old = old_ev2v[old_edge*2 + 0];
    LO const v1_old = old_ev2v[old_edge*2 + 1];

    LO const mid_vert = keys2midverts[key];
    LO const start = keys2prods[key];
    LO const end = keys2prods[key + 1] - 1;
    LO const new_e0 = prods2new[start];
    LO const new_e1 = prods2new[start+1];

    LO const v1_new_e0 = new_ev2v[new_e0*2 + 1];
    LO const v0_new_e1 = new_ev2v[new_e1*2 + 0];
    OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));

    //ctrl pts for e0
    {
      Real new_xi_start = 0.0;
      Real const cx0 = old_vertCtrlPts[v0_old*dim + 0];
      Real const cy0 = old_vertCtrlPts[v0_old*dim + 1];
      Real const cz0 = old_vertCtrlPts[v0_old*dim + 2];
      Real const cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
      Real const cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
      Real const cz1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 2];
      Real const cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
        (n_edge_pts-1)*dim + 0];
      Real const cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
        (n_edge_pts-1)*dim + 1];
      Real const cz2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
        (n_edge_pts-1)*dim + 2];
      Real const cx3 = old_vertCtrlPts[v1_old*dim + 0];
      Real const cy3 = old_vertCtrlPts[v1_old*dim + 1];
      Real const cz3 = old_vertCtrlPts[v1_old*dim + 2];

      Real const new_xi_3 = new_xi_start + 0.5;
      Real const new_cx3 = cx0*Bi(order, 0, new_xi_3) + cx1*Bi(order, 1, new_xi_3) +
        cx2*Bi(order, 2, new_xi_3) + cx3*Bi(order, 3, new_xi_3);
      Real const new_cy3 = cy0*Bi(order, 0, new_xi_3) + cy1*Bi(order, 1, new_xi_3) +
        cy2*Bi(order, 2, new_xi_3) + cy3*Bi(order, 3, new_xi_3);
      Real const new_cz3 = cz0*Bi(order, 0, new_xi_3) + cz1*Bi(order, 1, new_xi_3) +
        cz2*Bi(order, 2, new_xi_3) + cz3*Bi(order, 3, new_xi_3);
        
      vert_ctrlPts[mid_vert*1*dim + 0] = new_cx3;
      vert_ctrlPts[mid_vert*1*dim + 1] = new_cy3;
      vert_ctrlPts[mid_vert*1*dim + 2] = new_cz3;

      Real const old_xi_2 = xi_2_cube();
      Real const new_xi_2 = new_xi_start + old_xi_2/2.0;
      Real const new_px2 = cx0*Bi(order, 0, new_xi_2) + cx1*Bi(order, 1, new_xi_2) +
        cx2*Bi(order, 2, new_xi_2) + cx3*Bi(order, 3, new_xi_2);
      Real const new_py2 = cy0*Bi(order, 0, new_xi_2) + cy1*Bi(order, 1, new_xi_2) +
        cy2*Bi(order, 2, new_xi_2) + cy3*Bi(order, 3, new_xi_2);
      Real const new_pz2 = cz0*Bi(order, 0, new_xi_2) + cz1*Bi(order, 1, new_xi_2) +
        cz2*Bi(order, 2, new_xi_2) + cz3*Bi(order, 3, new_xi_2);

      Real const old_xi_1 = xi_1_cube();
      Real const new_xi_1 = new_xi_start + old_xi_1/2.0;
      Real const new_px1 = cx0*Bi(order, 0, new_xi_1) + cx1*Bi(order, 1, new_xi_1) +
        cx2*Bi(order, 2, new_xi_1) + cx3*Bi(order, 3, new_xi_1);
      Real const new_py1 = cy0*Bi(order, 0, new_xi_1) + cy1*Bi(order, 1, new_xi_1) +
        cy2*Bi(order, 2, new_xi_1) + cy3*Bi(order, 3, new_xi_1);
      Real const new_pz1 = cz0*Bi(order, 0, new_xi_1) + cz1*Bi(order, 1, new_xi_1) +
        cz2*Bi(order, 2, new_xi_1) + cz3*Bi(order, 3, new_xi_1);

      auto fx = vector_2(new_px1, new_px2);
      auto fy = vector_2(new_py1, new_py2);
      auto fz = vector_2(new_pz1, new_pz2);
      auto M1_inv = matrix_2x2(Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
          Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2));
      auto M2 = matrix_2x2(Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
          Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2));
      auto M1 = invert(M1_inv);
      auto cx = vector_2(cx0, new_cx3);
      auto cy = vector_2(cy0, new_cy3);
      auto cz = vector_2(cz0, new_cz3);

      auto Cx = M1*fx - M1*M2*cx;
      auto Cy = M1*fy - M1*M2*cy;
      auto Cz = M1*fz - M1*M2*cz;
      edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx[0];
      edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy[0];
      edge_ctrlPts[new_e0*n_edge_pts*dim + 2] = Cz[0];
      edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
      edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
      edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz[1];
    }

    //ctrl pts for e1
    {
      Real new_xi_start = 0.5;
      Real cx0 = old_vertCtrlPts[v0_old*dim + 0];
      Real cy0 = old_vertCtrlPts[v0_old*dim + 1];
      Real cz0 = old_vertCtrlPts[v0_old*dim + 2];
      Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
      Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
      Real cz1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 2];
      Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0];
      Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1];
      Real cz2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 2];
      Real cx3 = old_vertCtrlPts[v1_old*dim + 0];
      Real cy3 = old_vertCtrlPts[v1_old*dim + 1];
      Real cz3 = old_vertCtrlPts[v1_old*dim + 2];

      Real old_xi_2 = xi_2_cube();
      Real new_xi_2 = new_xi_start + old_xi_2/2.0;
      Real const new_px2 = cx0*Bi(order, 0, new_xi_2) + cx1*Bi(order, 1, new_xi_2) +
        cx2*Bi(order, 2, new_xi_2) + cx3*Bi(order, 3, new_xi_2);
      Real const new_py2 = cy0*Bi(order, 0, new_xi_2) + cy1*Bi(order, 1, new_xi_2) +
        cy2*Bi(order, 2, new_xi_2) + cy3*Bi(order, 3, new_xi_2);
      Real const new_pz2 = cz0*Bi(order, 0, new_xi_2) + cz1*Bi(order, 1, new_xi_2) +
        cz2*Bi(order, 2, new_xi_2) + cz3*Bi(order, 3, new_xi_2);

      Real old_xi_1 = xi_1_cube();
      Real new_xi_1 = new_xi_start + old_xi_1/2.0;
      Real const new_px1 = cx0*Bi(order, 0, new_xi_1) + cx1*Bi(order, 1, new_xi_1) +
        cx2*Bi(order, 2, new_xi_1) + cx3*Bi(order, 3, new_xi_1);
      Real const new_py1 = cy0*Bi(order, 0, new_xi_1) + cy1*Bi(order, 1, new_xi_1) +
        cy2*Bi(order, 2, new_xi_1) + cy3*Bi(order, 3, new_xi_1);
      Real const new_pz1 = cz0*Bi(order, 0, new_xi_1) + cz1*Bi(order, 1, new_xi_1) +
        cz2*Bi(order, 2, new_xi_1) + cz3*Bi(order, 3, new_xi_1);

      auto fx = vector_2(new_px1, new_px2);
      auto fy = vector_2(new_py1, new_py2);
      auto fz = vector_2(new_pz1, new_pz2);
      auto M1_inv = matrix_2x2(Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
          Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2));
      auto M2 = matrix_2x2(Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
          Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2));
      auto M1 = invert(M1_inv);
      auto cx = vector_2(vert_ctrlPts[mid_vert*1*dim + 0], cx3);
      auto cy = vector_2(vert_ctrlPts[mid_vert*1*dim + 1], cy3);
      auto cz = vector_2(vert_ctrlPts[mid_vert*1*dim + 2], cz3);

      auto Cx = M1*fx - M1*M2*cx;
      auto Cy = M1*fy - M1*M2*cy;
      auto Cz = M1*fz - M1*M2*cz;
      edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx[0];
      edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy[0];
      edge_ctrlPts[new_e1*n_edge_pts*dim + 2] = Cz[0];
      edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
      edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
      edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz[1];
    }
    //ctrl pts for edges on adjacent faces
    for (LO i = 0; i <= (end-start - 2); ++i) {
      LO const new_e2 = prods2new[start+2 + i];

      LO const v0_new_e2 = new_ev2v[new_e2*2 + 0];
      LO const v1_new_e2 = new_ev2v[new_e2*2 + 1];
      OMEGA_H_CHECK(v1_new_e2 == mid_vert);

      auto ab2b = new_verts2old_verts.ab2b;
      auto a2ab = new_verts2old_verts.a2ab;
      OMEGA_H_CHECK((a2ab[v0_new_e2+1] - a2ab[v0_new_e2]) == 1);
      LO old_vert_noKey = ab2b[a2ab[v0_new_e2]];

      LO old_face = -1;
      for (LO index = old_e2ef[old_edge]; index < old_e2ef[old_edge + 1];
          ++index) {
        LO const adj_face = old_ef2f[index];
        for (LO vert = 0; vert < 3; ++vert) {
          LO const vert_old_face = old_fv2v[adj_face*3 + vert];
          if (vert_old_face == old_vert_noKey) {
            old_face = adj_face;
            break;
          }
        }
        if (old_face > 0) {
          break;
        }
      }

      keys2old_faces_w[max_degree_key2oldface*key + i] = old_face;
      {
        auto const v0_old_face = old_fv2v[old_face*3];
        auto const v1 = old_fv2v[old_face*3 + 1];
        auto const v2 = old_fv2v[old_face*3 + 2];
        auto const old_face_e0 = old_fe2e[old_face*3 + 0];
        auto const old_face_e1 = old_fe2e[old_face*3 + 1];
        auto const old_face_e2 = old_fe2e[old_face*3 + 2];

        auto nodePts = cubic_noKeyEdge_xi_values(old_vert_noKey, v0_old_face, v1, v2,
            old_edge, old_face_e0, old_face_e1,
            old_face_e2);
        //get the interp points
        auto p1 = face_parametricToParent_3d(order, old_face, old_ev2v, old_fe2e,
            old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[0], nodePts[1], old_fv2v);
        auto p2 = face_parametricToParent_3d(order, old_face, old_ev2v, old_fe2e,
            old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[2], nodePts[3], old_fv2v);

        //use these as interp pts to find ctrl pts for the new mesh
        {
          Real cx0 = old_vertCtrlPts[old_vert_noKey*dim + 0];
          Real cy0 = old_vertCtrlPts[old_vert_noKey*dim + 1];
          Real cz0 = old_vertCtrlPts[old_vert_noKey*dim + 2];

          auto fx = vector_2(p1[0], p2[0]);
          auto fy = vector_2(p1[1], p2[1]);
          auto fz = vector_2(p1[2], p2[2]);
          Real old_xi_1 = xi_1_cube();
          Real old_xi_2 = xi_2_cube();
          auto M1_inv = matrix_2x2(Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
              Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2));
          auto M2 = matrix_2x2(Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
              Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2));
          auto M1 = invert(M1_inv);
          auto cx = vector_2(cx0, vert_ctrlPts[mid_vert*dim + 0]);
          auto cy = vector_2(cy0, vert_ctrlPts[mid_vert*dim + 1]);
          auto cz = vector_2(cz0, vert_ctrlPts[mid_vert*dim + 2]);

          auto Cx = M1*fx - M1*M2*cx;
          auto Cy = M1*fy - M1*M2*cy;
          auto Cz = M1*fz - M1*M2*cz;
          edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx[0];
          edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy[0];
          edge_ctrlPts[new_e2*n_edge_pts*dim + 2] = Cz[0];
          edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
          edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
          edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz[1];
        }
      }
    }
  };
  parallel_for(nkeys, std::move(create_crv_prod_edges));

  auto create_crv_same_edges = OMEGA_H_LAMBDA (LO old_edge) {
    if (old2new[old_edge] != -1) {
      LO new_edge = old2new[old_edge];
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*n_edge_pts*dim + d] =
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + d];
        edge_ctrlPts[new_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d] =
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d];
      }
    }
  };
  parallel_for(nold_edge, std::move(create_crv_same_edges));

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim);
  new_mesh->add_tag<Real>(0, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));

  //copy ctrl pts for same verts
  auto copy_sameCtrlPts = OMEGA_H_LAMBDA(LO i) {
    if (old_verts2new_verts[i] != -1) {
      LO new_vert = old_verts2new_verts[i];
      for (I8 d = 0; d < dim; ++d) {
        vert_ctrlPts[new_vert*dim + d] = old_vertCtrlPts[i*dim + d];
      }
    }
  };
  parallel_for(nold_verts, std::move(copy_sameCtrlPts), "copy same vtx ctrlPts");
  new_mesh->set_tag_for_ctrlPts(0, Reals(vert_ctrlPts));

  return LOs(keys2old_faces_w);
}

void create_curved_faces_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                           LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                           LOs old_verts2new_verts) {
  auto const order = mesh->get_max_order();
  OMEGA_H_CHECK(order >= 2);
  OMEGA_H_TIME_FUNCTION;
  auto const nold_face = old2new.size();
  auto const dim = mesh->dim();

  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_er2r = mesh->ask_up(1, 3).ab2b;
  auto const old_e2er = mesh->ask_up(1, 3).a2ab;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto const old_rv2v = mesh->ask_down(3, 0).ab2b;
  auto const old_re2e = mesh->ask_down(3, 1).ab2b;
  auto const old_rf2f = mesh->get_adj(3, 2).ab2b;
  auto const old_coords = mesh->coords();

  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_face_pts = mesh->n_internal_ctrlPts(2);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_fe2e = new_mesh->get_adj(2, 1).ab2b;
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_face = new_mesh->nfaces();
  auto const nnew_verts = new_mesh->nverts();
  auto const new_edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const new_vertCtrlPts = new_mesh->get_ctrlPts(0);
  OMEGA_H_CHECK(nnew_verts == (new_vertCtrlPts.size()/dim));

  Write<Real> face_ctrlPts(nnew_face*n_face_pts*dim);
  OMEGA_H_CHECK(order == 3);

  auto const nkeys = keys2edges.size();
  LO const max_degree_key2oldface = keys2old_faces.size()/nkeys;

  auto keys2nold_faces_w = Write<LO>(nkeys, -1);
  auto count_nold_faces = OMEGA_H_LAMBDA (LO key) {
    for (LO i = 0; i < max_degree_key2oldface; ++i) {
      if (keys2old_faces[max_degree_key2oldface*key + i] != -1) {
        atomic_increment(&keys2nold_faces_w[key]);
      }
    }
  };
  parallel_for(nkeys, count_nold_faces);
  auto keys2nold_faces = LOs(keys2nold_faces_w);

  auto create_crv_prod_faces = OMEGA_H_LAMBDA (LO key) {

    auto start = keys2prods[key];
    auto end = keys2prods[key + 1] - 1;

    LO const new_faces_per_old_face = 2;

    for (LO i = 0; i <= keys2nold_faces[key]; ++i) {
      auto new_f0 = prods2new[start + (i*new_faces_per_old_face) + 0];
      auto new_f1 = prods2new[start + (i*new_faces_per_old_face) + 1];

      auto old_face = keys2old_faces[max_degree_key2oldface*key + i];
      OMEGA_H_CHECK(old_face != -1);
      auto old_key_edge = keys2edges[key];

      LO const v0_old_face = old_fv2v[old_face*3 + 0];
      LO const v1_old_face = old_fv2v[old_face*3 + 1];
      LO const v2_old_face = old_fv2v[old_face*3 + 2];
      auto const old_face_e0 = old_fe2e[old_face*3];
      auto const old_face_e1 = old_fe2e[old_face*3 + 1];
      auto const old_face_e2 = old_fe2e[old_face*3 + 2];
      auto old_key_edge_v0 = old_ev2v[old_key_edge*2 + 0];
      auto old_key_edge_v1 = old_ev2v[old_key_edge*2 + 1];
      LO old_vert_noKey = -1;
      for (LO k = 0; k < 3; ++k) {
        auto old_face_vert = old_fv2v[old_face*3 + k];
        if ((old_face_vert != old_key_edge_v0) &&
            (old_face_vert != old_key_edge_v1)) {
          old_vert_noKey = old_face_vert;
          break;
        }
      }

      auto nodePts = cubic_faceSplittingEdge_xi_values
        (old_vert_noKey, v0_old_face, v1_old_face, v2_old_face, old_key_edge,
         old_face_e0, old_face_e1, old_face_e2, new_fv2v[new_f0*3 + 0], 
         new_fv2v[new_f0*3 + 1], new_fv2v[new_f0*3 + 2], new_fv2v[new_f1*3 + 0], 
         new_fv2v[new_f1*3 + 1], new_fv2v[new_f1*3 + 2],
         old_verts2new_verts[v0_old_face],
         old_verts2new_verts[v1_old_face], old_verts2new_verts[v2_old_face]);

      //for f0
      {
        //get the interp points
        auto p11 = face_parametricToParent_3d(order, old_face, old_ev2v, old_fe2e,
          old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[0], nodePts[1], old_fv2v);
        //use these as interp pts to find ctrl pt in new face
        auto newface_c11 = face_interpToCtrlPt_3d(order, new_f0, new_ev2v, new_fe2e,
          new_vertCtrlPts, new_edgeCtrlPts, p11, new_fv2v);
        for (LO k = 0; k < 3; ++k) {
          face_ctrlPts[new_f0*n_face_pts*dim + k] = newface_c11[k];
          //TODO order > 3; (multiple pts. per face)
        }

      }

      //for f1
      {
        auto p11 = face_parametricToParent_3d(order, old_face, old_ev2v, old_fe2e,
          old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[2], nodePts[3], old_fv2v);
        auto newface_c11 = face_interpToCtrlPt_3d(order, new_f1, new_ev2v, new_fe2e,
          new_vertCtrlPts, new_edgeCtrlPts, p11, new_fv2v);
        for (LO k = 0; k < 3; ++k) {
          face_ctrlPts[new_f1*n_face_pts*dim + k] = newface_c11[k];
        }

      }
    }
    for (LO i = (keys2nold_faces[key]+1)*new_faces_per_old_face; i <= (end-start); ++i) {
      auto old_key_edge = keys2edges[key];
      auto newface = prods2new[start + i];
      LO const newface_v0 = new_fv2v[newface*3 + 0];
      LO const newface_v1 = new_fv2v[newface*3 + 1];
      LO const newface_v2 = new_fv2v[newface*3 + 2];

      LO old_rgn = -1;
      {
        auto adjrgn_start = old_e2er[old_key_edge];
        auto adjrgn_end = old_e2er[old_key_edge + 1];
        for (LO count_adj_rgn = adjrgn_start; count_adj_rgn < adjrgn_end;
             ++count_adj_rgn) {
          auto adj_rgn = old_er2r[count_adj_rgn];
          auto adj_rgn_v0_old = old_rv2v[adj_rgn*4 + 0];
          auto adj_rgn_v1_old = old_rv2v[adj_rgn*4 + 1];
          auto adj_rgn_v2_old = old_rv2v[adj_rgn*4 + 2];
          auto adj_rgn_v3_old = old_rv2v[adj_rgn*4 + 3];
          auto adj_rgn_v0_new = old_verts2new_verts[adj_rgn_v0_old];
          auto adj_rgn_v1_new = old_verts2new_verts[adj_rgn_v1_old];
          auto adj_rgn_v2_new = old_verts2new_verts[adj_rgn_v2_old];
          auto adj_rgn_v3_new = old_verts2new_verts[adj_rgn_v3_old];

          LO count_matchverts2 = 0;
          if ((newface_v0 == adj_rgn_v0_new) || (newface_v0 == adj_rgn_v1_new)||
              (newface_v0 == adj_rgn_v2_new) || (newface_v0 == adj_rgn_v3_new)) {
            ++count_matchverts2;
          }
          if ((newface_v1 == adj_rgn_v0_new) || (newface_v1 == adj_rgn_v1_new)||
              (newface_v1 == adj_rgn_v2_new) || (newface_v1 == adj_rgn_v3_new)) {
            ++count_matchverts2;
          }
          if ((newface_v2 == adj_rgn_v0_new) || (newface_v2 == adj_rgn_v1_new)||
              (newface_v2 == adj_rgn_v2_new) || (newface_v2 == adj_rgn_v3_new)) {
            ++count_matchverts2;
          }
          if (count_matchverts2 == 2) {
            old_rgn = adj_rgn;
            break;
          }
        }
      }

      auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
      auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
      auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
      auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
      auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
      auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];
      auto nodePt = cubic_region_xi_values(
        old_key_edge, old_rgn_e0, old_rgn_e1, old_rgn_e2, old_rgn_e3,
        old_rgn_e4, old_rgn_e5);

      auto p11 = rgn_parametricToParent_3d(order, old_rgn, old_ev2v, old_rv2v,
          old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePt, old_re2e,
          old_rf2f);
      auto c11 = face_interpToCtrlPt_3d(order, newface, new_ev2v, new_fe2e,
          new_vertCtrlPts, new_edgeCtrlPts, p11, new_fv2v);
      for (LO j = 0; j < dim; ++j) {
        face_ctrlPts[newface*n_face_pts*dim + j] = c11[j];
      }
    }
  };
  parallel_for(nkeys, std::move(create_crv_prod_faces),
               "create_crv_prod_faces");

  auto create_crv_same_faces = OMEGA_H_LAMBDA (LO old_face) {
    if (old2new[old_face] != -1) {
      LO new_face = old2new[old_face];
      for (I8 d = 0; d < dim; ++d) {
        face_ctrlPts[new_face*n_face_pts*dim + d] = 
          old_faceCtrlPts[old_face*n_face_pts*dim + d];
      }
    }
  };
  parallel_for(nold_face, std::move(create_crv_same_faces),
               "create_crv_same_faces");

  new_mesh->add_tag<Real>(2, "bezier_pts", n_face_pts*dim);
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
