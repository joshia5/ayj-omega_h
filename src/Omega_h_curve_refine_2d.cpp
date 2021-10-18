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

LOs create_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges) {
  printf("in create curved edges fn\n");
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
  auto const order = mesh->get_max_order();
  auto const dim = mesh->dim();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<Real> vert_ctrlPts(nnew_verts*1*dim, INT8_MAX);

  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  LO max_degree_key2oldface = -1;
  {
    auto keyedges_noldfaces_w = Write<LO>(nold_edge, -1);
    auto count_oldfaces = OMEGA_H_LAMBDA (LO old_edge) {
      if (old2new[old_edge] == -1) {//old edge is key
        //get num up adj faces of key edge
        LO const num_adj_faces = old_e2ef[old_edge + 1] - old_e2ef[old_edge];
        keyedges_noldfaces_w[old_edge] = num_adj_faces;
      }
    };
    parallel_for(nold_edge, std::move(count_oldfaces));
    max_degree_key2oldface = get_max(LOs(keyedges_noldfaces_w));
  }
  printf("max key2oldface degree = %d\n",max_degree_key2oldface);

  Write<LO> count_key(1, 0);
  auto nkeys = keys2midverts.size();
  OMEGA_H_CHECK(order == 3);
  auto keys2old_faces_w = Write<LO>(max_degree_key2oldface*nkeys, -1);

  auto create_crv_edges = OMEGA_H_LAMBDA (LO old_edge) {
    LO const v0_old = old_ev2v[old_edge*2 + 0];
    LO const v1_old = old_ev2v[old_edge*2 + 1];

    if (old2new[old_edge] != -1) {
      LO new_edge = old2new[old_edge];
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*n_edge_pts*dim + d] = 
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + d];
        edge_ctrlPts[new_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d] =
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d];
      }
    }
    else {
      LO const key_id = count_key[0];
      LO const mid_vert = keys2midverts[key_id];
      LO const start = keys2prods[key_id];
      LO const end = keys2prods[key_id + 1] - 1;
      OMEGA_H_CHECK(old_edge == keys2edges[key_id]);
      if ((end-start) != 2) {
        printf("for old edge %d, new edges=%d\n", old_edge, end-start+1);
        OMEGA_H_CHECK((end-start) == 3);
      }
      LO const new_e0 = prods2new[start];
      LO const new_e1 = prods2new[start+1];

      LO const v1_new_e0 = new_ev2v[new_e0*2 + 1];
      LO const v0_new_e1 = new_ev2v[new_e1*2 + 0];
      OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));

      //ctrl pts for e0
      {
        //TODO order=3
        Real const cx0 = old_vertCtrlPts[v0_old*dim + 0];
        Real const cy0 = old_vertCtrlPts[v0_old*dim + 1];
        Real const cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real const cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real const cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 0];
        Real const cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 1];
        Real const cx3 = old_vertCtrlPts[v1_old*dim + 0];
        Real const cy3 = old_vertCtrlPts[v1_old*dim + 1];

        Real const new_xi_3 = 0.5;
        Real const new_cx3 = cx0*Bi(order, 0, new_xi_3) + cx1*Bi(order, 1, new_xi_3) +
                       cx2*Bi(order, 2, new_xi_3) + cx3*Bi(order, 3, new_xi_3);
        Real const new_cy3 = cy0*Bi(order, 0, new_xi_3) + cy1*Bi(order, 1, new_xi_3) +
                       cy2*Bi(order, 2, new_xi_3) + cy3*Bi(order, 3, new_xi_3);
        printf(" new mid interp pt %f, %f \n", new_cx3, new_cy3);
        vert_ctrlPts[mid_vert*1*dim + 0] = new_cx3;
        vert_ctrlPts[mid_vert*1*dim + 1] = new_cy3;

        Real const old_xi_2 = xi_2_cube();
        Real const new_xi_2 = old_xi_2/2.0;
        Real const new_px2 = cx0*Bi(order, 0, new_xi_2) + cx1*Bi(order, 1, new_xi_2) +
                       cx2*Bi(order, 2, new_xi_2) + cx3*Bi(order, 3, new_xi_2);
        Real const new_py2 = cy0*Bi(order, 0, new_xi_2) + cy1*Bi(order, 1, new_xi_2) +
                       cy2*Bi(order, 2, new_xi_2) + cy3*Bi(order, 3, new_xi_2);
        printf("edge 0 p2 %f, %f , xi2 %f\n", new_px2, new_py2, new_xi_2);

        Real const old_xi_1 = xi_1_cube();
        Real const new_xi_1 = old_xi_1/2.0;
        Real const new_px1 = cx0*Bi(order, 0, new_xi_1) + cx1*Bi(order, 1, new_xi_1) +
                       cx2*Bi(order, 2, new_xi_1) + cx3*Bi(order, 3, new_xi_1);
        Real const new_py1 = cy0*Bi(order, 0, new_xi_1) + cy1*Bi(order, 1, new_xi_1) +
                       cy2*Bi(order, 2, new_xi_1) + cy3*Bi(order, 3, new_xi_1);
        printf("edge 0 p1 %f, %f \n", new_px1, new_py1);

        auto fx = vector_2(new_px1, new_px2);
        auto fy = vector_2(new_py1, new_py2);
        auto M1_inv = matrix_2x2(Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
                            Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2));
        auto M2 = matrix_2x2(Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
                        Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2));
        auto M1 = invert(M1_inv);
        auto cx = vector_2(cx0, new_cx3);
        auto cy = vector_2(cy0, new_cy3);
        
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx[0];
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy[0];
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
        /*
        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});
        Matrix<2,2> M1_inv({Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
                            Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2)});
        Matrix<2,2> M2({Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
                        Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2)});
        auto M1 = invert(M1_inv);
        Matrix<2,1> cx({cx0, new_cx3});
        Matrix<2,1> cy({cy0, new_cy3});
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        printf("\n edge 0 c1 %f, %f \n", Cx(0,0), Cy(0,0));
        printf(" edge 0 c2 %f, %f\n", Cx(1,0), Cy(1,0));
        */

      }
      //ctrl pts for e1
      {
        //for 2d mesh for now, order=3
        Real cx0 = old_vertCtrlPts[v0_old*dim + 0];
        Real cy0 = old_vertCtrlPts[v0_old*dim + 1];
        Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0];
        Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1];
        Real cx3 = old_vertCtrlPts[v1_old*dim + 0];
        Real cy3 = old_vertCtrlPts[v1_old*dim + 1];

        Real old_xi_2 = xi_2_cube();
        Real new_xi_2 = 0.5 + old_xi_2/2.0;
        Real const new_px2 = cx0*Bi(order, 0, new_xi_2) + cx1*Bi(order, 1, new_xi_2) +
                       cx2*Bi(order, 2, new_xi_2) + cx3*Bi(order, 3, new_xi_2);
        Real const new_py2 = cy0*Bi(order, 0, new_xi_2) + cy1*Bi(order, 1, new_xi_2) +
                       cy2*Bi(order, 2, new_xi_2) + cy3*Bi(order, 3, new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = 0.5 + old_xi_1/2.0;
        Real const new_px1 = cx0*Bi(order, 0, new_xi_1) + cx1*Bi(order, 1, new_xi_1) +
                       cx2*Bi(order, 2, new_xi_1) + cx3*Bi(order, 3, new_xi_1);
        Real const new_py1 = cy0*Bi(order, 0, new_xi_1) + cy1*Bi(order, 1, new_xi_1) +
                       cy2*Bi(order, 2, new_xi_1) + cy3*Bi(order, 3, new_xi_1);

        auto cx = vector_2(vert_ctrlPts[mid_vert*1*dim + 0], cx3);
        auto cy = vector_2(vert_ctrlPts[mid_vert*1*dim + 1], cy3);
        auto fx = vector_2(new_px1, new_px2);
        auto fy = vector_2(new_py1, new_py2);
        auto M1_inv = matrix_2x2(Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
                            Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2));
        auto M2 = matrix_2x2(Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
                        Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2));
        auto M1 = invert(M1_inv);
        
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx[0];
        edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy[0];
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
        /*
        Matrix<2,1> cx({vert_ctrlPts[mid_vert*1*dim + 0], cx3});
        Matrix<2,1> cy({vert_ctrlPts[mid_vert*1*dim + 1], cy3});
                        
        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({Bi(order, 1, old_xi_1), Bi(order, 2, old_xi_1),
                            Bi(order, 1, old_xi_2), Bi(order, 2, old_xi_2)});
        Matrix<2,2> M2({Bi(order, 0, old_xi_1), Bi(order, 3, old_xi_1),
                        Bi(order, 0, old_xi_2), Bi(order, 3, old_xi_2)});

        auto M1 = invert(M1_inv);
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        */
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
          LO adj_face = old_ef2f[index];
          for (LO vert = 0; vert < 3; ++vert) {
            LO vert_old_face = old_fv2v[adj_face*3 + vert];
            if (vert_old_face == old_vert_noKey) {
              old_face = adj_face;
              break;
            }
          }
          if (old_face > 0) {
            break;
          }
        }

        printf("For oldedge %d with key id %d, found oldface %d\n", old_edge, key_id, old_face);
        keys2old_faces_w[max_degree_key2oldface*key_id + i] = old_face;
        {
          auto const v0_old_face = old_fv2v[old_face*3];
          auto const v1 = old_fv2v[old_face*3 + 1];
          auto const v2 = old_fv2v[old_face*3 + 2];
          auto const old_face_e0 = old_fe2e[old_face*3];
          auto const old_face_e1 = old_fe2e[old_face*3 + 1];
          auto const old_face_e2 = old_fe2e[old_face*3 + 2];

          auto nodePts = cubic_noKeyEdge_xi_values(old_vert_noKey, v0_old_face, v1, v2,
                                            old_edge, old_face_e0, old_face_e1,
                                            old_face_e2);

          auto p1 = face_parametricToParent_2d(order, old_face, old_ev2v, old_fe2e,
              old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[0], nodePts[1], old_fv2v);
          auto p2 = face_parametricToParent_2d(order, old_face, old_ev2v, old_fe2e,
              old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[2], nodePts[3], old_fv2v);
          //use these as interp pts to find ctrl pts in new mesh
          {
            Real cx0 = old_vertCtrlPts[old_vert_noKey*dim + 0];
            Real cy0 = old_vertCtrlPts[old_vert_noKey*dim + 1];

            auto cx = vector_2(cx0, vert_ctrlPts[mid_vert*1*dim + 0]);
            auto cy = vector_2(cy0, vert_ctrlPts[mid_vert*1*dim + 1]);
            auto fx = vector_2(p1[0], p2[0]);
            auto fy = vector_2(p1[1], p2[1]);
            auto M1_inv = matrix_2x2(Bi(order, 1, xi_1_cube()), Bi(order, 2, xi_1_cube()),
                Bi(order, 1, xi_2_cube()), Bi(order, 2, xi_2_cube()));
            auto M2 = matrix_2x2(Bi(order, 0, xi_1_cube()), Bi(order, 3, xi_1_cube()),
                Bi(order, 0, xi_2_cube()), Bi(order, 3, xi_2_cube()));
            auto M1 = invert(M1_inv);

            auto Cx = M1*fx - M1*M2*cx;
            auto Cy = M1*fy - M1*M2*cy;
            edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx[0];
            edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy[0];
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx[1];
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy[1];
        /*
            Matrix<2,1> cx({cx0, vert_ctrlPts[mid_vert*1*dim + 0]});
            Matrix<2,1> cy({cy0, vert_ctrlPts[mid_vert*1*dim + 1]});

            Matrix<2,1> fx({p1[0], p2[0]});
            Matrix<2,1> fy({p1[1], p2[1]});

            Matrix<2,2> M1_inv({Bi(order, 1, xi_1_cube()), Bi(order, 2, xi_1_cube()),
                Bi(order, 1, xi_2_cube()), Bi(order, 2, xi_2_cube())});
            Matrix<2,2> M2({Bi(order, 0,  xi_1_cube()), Bi(order, 3,  xi_1_cube()),
                Bi(order, 0,  xi_2_cube()), Bi(order, 3,  xi_2_cube())});

            auto M1 = invert(M1_inv);
            auto Cx = M1*fx - M1*M2*cx;
            auto Cy = M1*fy - M1*M2*cy;
            edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
            */
          }
        }
      }
      atomic_fetch_add(&count_key[0], 1);
    }
  };
  parallel_for(nold_edge, std::move(create_crv_edges));

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

void create_curved_faces_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts) {
  printf("in create curved faces fn\n");
  OMEGA_H_TIME_FUNCTION;
  auto const nold_face = old2new.size();
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const order = mesh->get_max_order();
  auto const dim = mesh->dim();
  //auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
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

  Write<Real> face_ctrlPts(nnew_face*n_face_pts*dim, INT8_MAX);
  OMEGA_H_CHECK(order == 3);

  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  auto const nkeys = keys2edges.size();
  LO const max_degree_key2oldface = keys2old_faces.size()/nkeys;

  auto create_crv_prod_faces = OMEGA_H_LAMBDA (LO key) {

    auto start = keys2prods[key];
    auto end = keys2prods[key + 1] - 1;
    if ((end-start) != 1) {
      printf("key %d split into %d faces\n", key, end-start+1);
      OMEGA_H_CHECK((end-start) == 3);
    }

    LO const new_faces_per_old_face = 2;
    OMEGA_H_CHECK (((end-start-1) % 2) == 0);

    for (LO i = 0; i <= (end-start-1)/new_faces_per_old_face; ++i) {
      auto new_f0 = prods2new[start + (i*new_faces_per_old_face) + 0];
      auto new_f1 = prods2new[start + (i*new_faces_per_old_face) + 1];

      auto old_face = keys2old_faces[max_degree_key2oldface*key + i];
      if (old_face != -1) {
        auto old_key_edge = keys2edges[key];
        printf("for oldedge %d, oldface is %d\n", old_key_edge, old_face);

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
        printf("for old face %d , oldKeyEdge %d, found old no-key vert %d\n",
            old_face, old_key_edge, old_vert_noKey);

        printf("ok1, new faces f0 f1 %d %d oldface %d\n", new_f0, new_f1, old_face);
        auto nodePts = cubic_faceSplittingEdge_xi_values
          (old_vert_noKey, v0_old_face, v1_old_face, v2_old_face, old_key_edge,
           old_face_e0, old_face_e1, old_face_e2, new_fv2v[new_f0*3 + 0], 
           new_fv2v[new_f0*3 + 1], new_fv2v[new_f0*3 + 2], new_fv2v[new_f1*3 + 0], 
           new_fv2v[new_f1*3 + 1], new_fv2v[new_f1*3 + 2],
           old_verts2new_verts[v0_old_face],
           old_verts2new_verts[v1_old_face], old_verts2new_verts[v2_old_face]);
        printf("ok2\n");

        //for f0
        {
          //get the interp point
          auto p11 = face_parametricToParent_2d(order, old_face, old_ev2v, old_fe2e,
            old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[0], nodePts[1], old_fv2v);
          auto newface_c11 = face_interpToCtrlPt_2d(order, new_f0, new_ev2v, new_fe2e,
            new_vertCtrlPts, new_edgeCtrlPts, p11, new_fv2v);
          for (LO k = 0; k < dim; ++k) {
            face_ctrlPts[new_f0*n_face_pts*dim + k] = newface_c11[k];
          }
        }

        //for f1
        {
          //get the interp point
          auto p11 = face_parametricToParent_2d(order, old_face, old_ev2v, old_fe2e,
            old_vertCtrlPts, old_edgeCtrlPts, old_faceCtrlPts, nodePts[2], nodePts[3], old_fv2v);
          auto newface_c11 = face_interpToCtrlPt_2d(order, new_f1, new_ev2v, new_fe2e,
            new_vertCtrlPts, new_edgeCtrlPts, p11, new_fv2v);
          for (LO k = 0; k < dim; ++k) {
            face_ctrlPts[new_f1*n_face_pts*dim + k] = newface_c11[k];
          }
        }
      }
    }
  };
  parallel_for(nkeys, std::move(create_crv_prod_faces),
               "create_crv_prod_faces");

  auto create_crv_same_faces = OMEGA_H_LAMBDA (LO old_face) {
    if (old2new[old_face] != -1) {
      LO new_face = old2new[old_face];
      printf("old face %d same as new face %d\n", old_face, new_face);
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
