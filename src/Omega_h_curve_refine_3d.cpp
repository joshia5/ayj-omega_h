#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"

namespace Omega_h {

LOs create_curved_verts_and_edges_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts) {
  printf("in 3d create curved edges fn\n");
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
  //TODO optimize

  Write<LO> count_key(1, 0);
  auto nkeys = keys2midverts.size();
  OMEGA_H_CHECK(order == 3);
  auto keys2old_faces_w = Write<LO>(2*nkeys, -1);

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
      if ((end-start) != 2) {
        printf("for old edge %d, new edges=%d\n", old_edge, end-start+1);
        OMEGA_H_CHECK((end-start) == 3);
      }
      //TODO generalize this if more than 4 new edges are produced
      LO const new_e0 = prods2new[start];
      LO const new_e1 = prods2new[start+1];
      LO const new_e2 = prods2new[start+2];

      LO const v1_new_e0 = new_ev2v[new_e0*2 + 1];
      LO const v0_new_e1 = new_ev2v[new_e1*2 + 0];
      OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));

      LO const v0_new_e2 = new_ev2v[new_e2*2 + 0];
      LO const v1_new_e2 = new_ev2v[new_e2*2 + 1];
      OMEGA_H_CHECK(v1_new_e2 == mid_vert);

      //ctrl pts for e0
      {
        //TODO use vectors
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

        Real const new_xi_3 = 0.5;
        Real const new_cx3 = cx0*B0_cube(new_xi_3) + cx1*B1_cube(new_xi_3) +
                       cx2*B2_cube(new_xi_3) + cx3*B3_cube(new_xi_3);
        Real const new_cy3 = cy0*B0_cube(new_xi_3) + cy1*B1_cube(new_xi_3) +
                       cy2*B2_cube(new_xi_3) + cy3*B3_cube(new_xi_3);
        Real const new_cz3 = cz0*B0_cube(new_xi_3) + cz1*B1_cube(new_xi_3) +
                       cz2*B2_cube(new_xi_3) + cz3*B3_cube(new_xi_3);
        printf(" new mid interp pt %f, %f \n", new_cx3, new_cy3);
        vert_ctrlPts[mid_vert*1*dim + 0] = new_cx3;
        vert_ctrlPts[mid_vert*1*dim + 1] = new_cy3;
        vert_ctrlPts[mid_vert*1*dim + 2] = new_cz3;

        Real const old_xi_2 = xi_2_cube();
        Real const new_xi_2 = old_xi_2/2.0;
        Real const new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real const new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);
        Real const new_pz2 = cz0*B0_cube(new_xi_2) + cz1*B1_cube(new_xi_2) +
                       cz2*B2_cube(new_xi_2) + cz3*B3_cube(new_xi_2);
        printf("edge 0 p2 %f, %f , xi2 %f\n", new_px2, new_py2, new_xi_2);

        Real const old_xi_1 = xi_1_cube();
        Real const new_xi_1 = old_xi_1/2.0;
        Real const new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real const new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);
        Real const new_pz1 = cz0*B0_cube(new_xi_1) + cz1*B1_cube(new_xi_1) +
                       cz2*B2_cube(new_xi_1) + cz3*B3_cube(new_xi_1);
        printf("edge 0 p1 %f, %f \n", new_px1, new_py1);

        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});
        Matrix<2,1> fz({new_pz1, new_pz2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        auto M1 = invert(M1_inv);

        Matrix<2,1> cx({cx0, new_cx3});
        Matrix<2,1> cy({cy0, new_cy3});
        Matrix<2,1> cz({cz0, new_cz3});
        
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        auto Cz = M1*fz - M1*M2*cz;
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 2] = Cz(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz(1,0);
        printf("\n edge 0 c1 %f, %f \n", Cx(0,0), Cy(0,0));
        printf(" edge 0 c2 %f, %f\n", Cx(1,0), Cy(1,0));

      }
      //ctrl pts for e1
      {
        //for 2d mesh for now, order=3
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
        Real new_xi_2 = 0.5 + old_xi_2/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);
        Real new_pz2 = cz0*B0_cube(new_xi_2) + cz1*B1_cube(new_xi_2) +
                       cz2*B2_cube(new_xi_2) + cz3*B3_cube(new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = 0.5 + old_xi_1/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);
        Real new_pz1 = cz0*B0_cube(new_xi_1) + cz1*B1_cube(new_xi_1) +
                       cz2*B2_cube(new_xi_1) + cz3*B3_cube(new_xi_1);

        Matrix<2,1> cx({vert_ctrlPts[mid_vert*1*dim + 0], cx3});
        Matrix<2,1> cy({vert_ctrlPts[mid_vert*1*dim + 1], cy3});
        Matrix<2,1> cz({vert_ctrlPts[mid_vert*1*dim + 2], cz3});
                        
        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});
        Matrix<2,1> fz({new_pz1, new_pz2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        auto M1 = invert(M1_inv);
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        auto Cz = M1*fz - M1*M2*cz;
        edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + 2] = Cz(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz(1,0);
      }
      //ctrl pts for e2
      {
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

        printf("For old edge %d, found old face %d\n", old_edge, old_face);
        keys2old_faces_w[2*key_id + 0] = old_face;
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

          I8 e0_flip = -1;
          I8 e1_flip = -1;
          I8 e2_flip = -1;
          auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
          auto e0v1 = old_ev2v[old_face_e0*2 + 1];
          auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
          auto e1v1 = old_ev2v[old_face_e1*2 + 1];
          auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
          auto e2v1 = old_ev2v[old_face_e2*2 + 1];
          if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
            e0_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
          }
          if ((e1v0_old_face == v2) && (e1v1 == v1)) {
            e1_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
          }
          if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
            e2_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
          }

          Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
          Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
          Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
          Real cx30 = old_vertCtrlPts[v1*dim + 0];
          Real cy30 = old_vertCtrlPts[v1*dim + 1];
          Real cz30 = old_vertCtrlPts[v1*dim + 2];
          Real cx03 = old_vertCtrlPts[v2*dim + 0];
          Real cy03 = old_vertCtrlPts[v2*dim + 1];
          Real cz03 = old_vertCtrlPts[v2*dim + 2];

          auto pts_per_edge = n_edge_pts;
          Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
          Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
          Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
          Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e0_flip > 0) {
            auto tempx = cx10;
            auto tempy = cy10;
            auto tempz = cz10;
            cx10 = cx20;
            cy10 = cy20;
            cz10 = cz20;
            cx20 = tempx;
            cy20 = tempy;
            cz20 = tempz;
          }

          Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
          Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
          Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
          Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e1_flip > 0) {
            auto tempx = cx21;
            auto tempy = cy21;
            auto tempz = cz21;
            cx21 = cx12;
            cy21 = cy12;
            cz21 = cz12;
            cx12 = tempx;
            cy12 = tempy;
            cz12 = tempz;
          }

          Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
          Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
          Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
          Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e2_flip > 0) {
            auto tempx = cx02;
            auto tempy = cy02;
            auto tempz = cz02;
            cx02 = cx01;
            cy02 = cy01;
            cz02 = cz01;
            cx01 = tempx;
            cy01 = tempy;
            cz01 = tempz;
          }

          Real cx11 = old_faceCtrlPts[old_face*dim + 0];
          Real cy11 = old_faceCtrlPts[old_face*dim + 1];
          Real cz11 = old_faceCtrlPts[old_face*dim + 2];

          //get the interp points
          auto p1_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
                      cx10*B10_cube(nodePts[0], nodePts[1]) +
                      cx20*B20_cube(nodePts[0], nodePts[1]) +
                      cx30*B30_cube(nodePts[0], nodePts[1]) +
                      cx21*B21_cube(nodePts[0], nodePts[1]) +
                      cx12*B12_cube(nodePts[0], nodePts[1]) +
                      cx03*B03_cube(nodePts[0], nodePts[1]) +
                      cx02*B02_cube(nodePts[0], nodePts[1]) +
                      cx01*B01_cube(nodePts[0], nodePts[1]) +
                      cx11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
                      cy10*B10_cube(nodePts[0], nodePts[1]) +
                      cy20*B20_cube(nodePts[0], nodePts[1]) +
                      cy30*B30_cube(nodePts[0], nodePts[1]) +
                      cy21*B21_cube(nodePts[0], nodePts[1]) +
                      cy12*B12_cube(nodePts[0], nodePts[1]) +
                      cy03*B03_cube(nodePts[0], nodePts[1]) +
                      cy02*B02_cube(nodePts[0], nodePts[1]) +
                      cy01*B01_cube(nodePts[0], nodePts[1]) +
                      cy11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_z = cz00*B00_cube(nodePts[0], nodePts[1]) +
                      cz10*B10_cube(nodePts[0], nodePts[1]) +
                      cz20*B20_cube(nodePts[0], nodePts[1]) +
                      cz30*B30_cube(nodePts[0], nodePts[1]) +
                      cz21*B21_cube(nodePts[0], nodePts[1]) +
                      cz12*B12_cube(nodePts[0], nodePts[1]) +
                      cz03*B03_cube(nodePts[0], nodePts[1]) +
                      cz02*B02_cube(nodePts[0], nodePts[1]) +
                      cz01*B01_cube(nodePts[0], nodePts[1]) +
                      cz11*B11_cube(nodePts[0], nodePts[1]);
          auto p2_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
                      cx10*B10_cube(nodePts[2], nodePts[3]) +
                      cx20*B20_cube(nodePts[2], nodePts[3]) +
                      cx30*B30_cube(nodePts[2], nodePts[3]) +
                      cx21*B21_cube(nodePts[2], nodePts[3]) +
                      cx12*B12_cube(nodePts[2], nodePts[3]) +
                      cx03*B03_cube(nodePts[2], nodePts[3]) +
                      cx02*B02_cube(nodePts[2], nodePts[3]) +
                      cx01*B01_cube(nodePts[2], nodePts[3]) +
                      cx11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
                      cy10*B10_cube(nodePts[2], nodePts[3]) +
                      cy20*B20_cube(nodePts[2], nodePts[3]) +
                      cy30*B30_cube(nodePts[2], nodePts[3]) +
                      cy21*B21_cube(nodePts[2], nodePts[3]) +
                      cy12*B12_cube(nodePts[2], nodePts[3]) +
                      cy03*B03_cube(nodePts[2], nodePts[3]) +
                      cy02*B02_cube(nodePts[2], nodePts[3]) +
                      cy01*B01_cube(nodePts[2], nodePts[3]) +
                      cy11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_z = cz00*B00_cube(nodePts[2], nodePts[3]) +
                      cz10*B10_cube(nodePts[2], nodePts[3]) +
                      cz20*B20_cube(nodePts[2], nodePts[3]) +
                      cz30*B30_cube(nodePts[2], nodePts[3]) +
                      cz21*B21_cube(nodePts[2], nodePts[3]) +
                      cz12*B12_cube(nodePts[2], nodePts[3]) +
                      cz03*B03_cube(nodePts[2], nodePts[3]) +
                      cz02*B02_cube(nodePts[2], nodePts[3]) +
                      cz01*B01_cube(nodePts[2], nodePts[3]) +
                      cz11*B11_cube(nodePts[2], nodePts[3]);

          printf("for e2 p1 is %f %f, p2 is %f %f\n", p1_x, p1_y, p2_x, p2_y);
          //use these as interp pts to find ctrl pts in new mesh
          {
            Real cx0 = old_vertCtrlPts[old_vert_noKey*dim + 0];
            Real cy0 = old_vertCtrlPts[old_vert_noKey*dim + 1];
            Real cz0 = old_vertCtrlPts[old_vert_noKey*dim + 2];

            Matrix<2,1> cx({cx0, vert_ctrlPts[mid_vert*1*dim + 0]});
            Matrix<2,1> cy({cy0, vert_ctrlPts[mid_vert*1*dim + 1]});
            Matrix<2,1> cz({cz0, vert_ctrlPts[mid_vert*1*dim + 2]});

            Matrix<2,1> fx({p1_x, p2_x});
            Matrix<2,1> fy({p1_y, p2_y});
            Matrix<2,1> fz({p1_z, p2_z});

            Matrix<2,2> M1_inv({B1_cube(xi_1_cube()), B2_cube(xi_1_cube()), B1_cube(xi_2_cube()),
                B2_cube(xi_2_cube())});
            Matrix<2,2> M2({B0_cube(xi_1_cube()), B3_cube(xi_1_cube()), B0_cube(xi_2_cube()),
                B3_cube(xi_2_cube())});

            auto M1 = invert(M1_inv);
            auto Cx = M1*fx - M1*M2*cx;
            auto Cy = M1*fy - M1*M2*cy;
            auto Cz = M1*fy - M1*M2*cz;
            edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + 2] = Cz(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz(1,0);
          }
        }
      }
      //ctrl pts for e3
      if (end-start == 3) {
        LO const new_e3 = prods2new[start+3];
        LO const v0_new_e3 = new_ev2v[new_e3*2 + 0];
        LO const v1_new_e3 = new_ev2v[new_e3*2 + 1];
        OMEGA_H_CHECK(v1_new_e3 == mid_vert);

        auto ab2b = new_verts2old_verts.ab2b;
        auto a2ab = new_verts2old_verts.a2ab;
        OMEGA_H_CHECK((a2ab[v0_new_e3+1] - a2ab[v0_new_e3]) == 1);
        LO old_vert_noKey = ab2b[a2ab[v0_new_e3]];

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

        printf("For old edge %d, found old face %d\n", old_edge, old_face);
        keys2old_faces_w[2*key_id + 1] = old_face;
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

          I8 e0_flip = -1;
          I8 e1_flip = -1;
          I8 e2_flip = -1;
          auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
          auto e0v1 = old_ev2v[old_face_e0*2 + 1];
          auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
          auto e1v1 = old_ev2v[old_face_e1*2 + 1];
          auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
          auto e2v1 = old_ev2v[old_face_e2*2 + 1];
          if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
            e0_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
          }
          if ((e1v0_old_face == v2) && (e1v1 == v1)) {
            e1_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
          }
          if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
            e2_flip = 1;
          }
          else {
            OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
          }

          Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
          Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
          Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
          Real cx30 = old_vertCtrlPts[v1*dim + 0];
          Real cy30 = old_vertCtrlPts[v1*dim + 1];
          Real cz30 = old_vertCtrlPts[v1*dim + 2];
          Real cx03 = old_vertCtrlPts[v2*dim + 0];
          Real cy03 = old_vertCtrlPts[v2*dim + 1];
          Real cz03 = old_vertCtrlPts[v2*dim + 2];

          auto pts_per_edge = n_edge_pts;
          Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
          Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
          Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
          Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e0_flip > 0) {
            swap2(cx10, cx20);
            swap2(cy10, cy20);
            swap2(cz10, cz20);
          }

          Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
          Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
          Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
          Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e1_flip > 0) {
            swap2(cx12, cx21);
            swap2(cy12, cy21);
            swap2(cz12, cz21);
          }

          Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
          Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
          Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
          Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e2_flip > 0) {
            swap2(cx02, cx01);
            swap2(cy02, cy01);
            swap2(cz02, cz01);
          }

          Real cx11 = old_faceCtrlPts[old_face*dim + 0];
          Real cy11 = old_faceCtrlPts[old_face*dim + 1];
          Real cz11 = old_faceCtrlPts[old_face*dim + 2];

          //get the interp points
          auto p1_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
                      cx10*B10_cube(nodePts[0], nodePts[1]) +
                      cx20*B20_cube(nodePts[0], nodePts[1]) +
                      cx30*B30_cube(nodePts[0], nodePts[1]) +
                      cx21*B21_cube(nodePts[0], nodePts[1]) +
                      cx12*B12_cube(nodePts[0], nodePts[1]) +
                      cx03*B03_cube(nodePts[0], nodePts[1]) +
                      cx02*B02_cube(nodePts[0], nodePts[1]) +
                      cx01*B01_cube(nodePts[0], nodePts[1]) +
                      cx11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
                      cy10*B10_cube(nodePts[0], nodePts[1]) +
                      cy20*B20_cube(nodePts[0], nodePts[1]) +
                      cy30*B30_cube(nodePts[0], nodePts[1]) +
                      cy21*B21_cube(nodePts[0], nodePts[1]) +
                      cy12*B12_cube(nodePts[0], nodePts[1]) +
                      cy03*B03_cube(nodePts[0], nodePts[1]) +
                      cy02*B02_cube(nodePts[0], nodePts[1]) +
                      cy01*B01_cube(nodePts[0], nodePts[1]) +
                      cy11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_z = cz00*B00_cube(nodePts[0], nodePts[1]) +
                      cz10*B10_cube(nodePts[0], nodePts[1]) +
                      cz20*B20_cube(nodePts[0], nodePts[1]) +
                      cz30*B30_cube(nodePts[0], nodePts[1]) +
                      cz21*B21_cube(nodePts[0], nodePts[1]) +
                      cz12*B12_cube(nodePts[0], nodePts[1]) +
                      cz03*B03_cube(nodePts[0], nodePts[1]) +
                      cz02*B02_cube(nodePts[0], nodePts[1]) +
                      cz01*B01_cube(nodePts[0], nodePts[1]) +
                      cz11*B11_cube(nodePts[0], nodePts[1]);
          auto p2_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
                      cx10*B10_cube(nodePts[2], nodePts[3]) +
                      cx20*B20_cube(nodePts[2], nodePts[3]) +
                      cx30*B30_cube(nodePts[2], nodePts[3]) +
                      cx21*B21_cube(nodePts[2], nodePts[3]) +
                      cx12*B12_cube(nodePts[2], nodePts[3]) +
                      cx03*B03_cube(nodePts[2], nodePts[3]) +
                      cx02*B02_cube(nodePts[2], nodePts[3]) +
                      cx01*B01_cube(nodePts[2], nodePts[3]) +
                      cx11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
                      cy10*B10_cube(nodePts[2], nodePts[3]) +
                      cy20*B20_cube(nodePts[2], nodePts[3]) +
                      cy30*B30_cube(nodePts[2], nodePts[3]) +
                      cy21*B21_cube(nodePts[2], nodePts[3]) +
                      cy12*B12_cube(nodePts[2], nodePts[3]) +
                      cy03*B03_cube(nodePts[2], nodePts[3]) +
                      cy02*B02_cube(nodePts[2], nodePts[3]) +
                      cy01*B01_cube(nodePts[2], nodePts[3]) +
                      cy11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_z = cz00*B00_cube(nodePts[2], nodePts[3]) +
                      cz10*B10_cube(nodePts[2], nodePts[3]) +
                      cz20*B20_cube(nodePts[2], nodePts[3]) +
                      cz30*B30_cube(nodePts[2], nodePts[3]) +
                      cz21*B21_cube(nodePts[2], nodePts[3]) +
                      cz12*B12_cube(nodePts[2], nodePts[3]) +
                      cz03*B03_cube(nodePts[2], nodePts[3]) +
                      cz02*B02_cube(nodePts[2], nodePts[3]) +
                      cz01*B01_cube(nodePts[2], nodePts[3]) +
                      cz11*B11_cube(nodePts[2], nodePts[3]);

          printf("for e3 p1 is %f %f, p2 is %f %f\n", p1_x, p1_y, p2_x, p2_y);
          //use these as interp pts to find ctrl pts in new mesh
          {
            Real cx0 = old_vertCtrlPts[old_vert_noKey*dim + 0];
            Real cy0 = old_vertCtrlPts[old_vert_noKey*dim + 1];
            Real cz0 = old_vertCtrlPts[old_vert_noKey*dim + 2];

            Matrix<2,1> cx({cx0, vert_ctrlPts[mid_vert*1*dim + 0]});
            Matrix<2,1> cy({cy0, vert_ctrlPts[mid_vert*1*dim + 1]});
            Matrix<2,1> cz({cz0, vert_ctrlPts[mid_vert*1*dim + 2]});

            Matrix<2,1> fx({p1_x, p2_x});
            Matrix<2,1> fy({p1_y, p2_y});
            Matrix<2,1> fz({p1_z, p2_z});

            Matrix<2,2> M1_inv({B1_cube(xi_1_cube()), B2_cube(xi_1_cube()), B1_cube(xi_2_cube()),
                B2_cube(xi_2_cube())});
            Matrix<2,2> M2({B0_cube(xi_1_cube()), B3_cube(xi_1_cube()), B0_cube(xi_2_cube()),
                B3_cube(xi_2_cube())});

            auto M1 = invert(M1_inv);
            auto Cx = M1*fx - M1*M2*cx;
            auto Cy = M1*fy - M1*M2*cy;
            auto Cz = M1*fy - M1*M2*cz;
            edge_ctrlPts[new_e3*n_edge_pts*dim + 0] = Cx(0,0);
            edge_ctrlPts[new_e3*n_edge_pts*dim + 1] = Cy(0,0);
            edge_ctrlPts[new_e3*n_edge_pts*dim + 2] = Cz(0,0);
            edge_ctrlPts[new_e3*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
            edge_ctrlPts[new_e3*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
            edge_ctrlPts[new_e3*n_edge_pts*dim + (n_edge_pts-1)*dim + 2] = Cz(1,0);
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
        OMEGA_H_CHECK(vert_ctrlPts[new_vert*dim + d] = INT8_MAX);
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
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto const n_face_pts = mesh->n_internal_ctrlPts(2);

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_fe2e = new_mesh->get_adj(2, 1).ab2b;
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_face = new_mesh->nfaces();
  //auto const nnew_edge = new_mesh->nedges();
  auto const nnew_verts = new_mesh->nverts();
  auto const new_edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const new_vertCtrlPts = new_mesh->get_ctrlPts(0);
  OMEGA_H_CHECK(nnew_verts == (new_vertCtrlPts.size()/dim));

  Write<Real> face_ctrlPts(nnew_face*n_face_pts*dim, INT8_MAX);
  OMEGA_H_CHECK(order == 3);

  //auto max_vert_old2new = get_max(old_verts2new_verts);
  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  auto const nkeys = keys2edges.size();
  auto create_crv_prod_faces = OMEGA_H_LAMBDA (LO key) {

    auto start = keys2prods[key];
    auto end = keys2prods[key + 1] - 1;
    if ((end-start) != 1) {
      printf("key %d split into %d faces\n", key, end-start+1);
      OMEGA_H_CHECK((end-start) == 3);
    }
    //TODO case where 1 key makes more than 4 faces
    {
      auto new_f0 = prods2new[start];
      auto new_f1 = prods2new[start + 1];

      auto old_face = keys2old_faces[2*key + 0];
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
      printf("for old face %d , oldKeyEdge %d, found old no-key vert %d\n",
          old_face, old_key_edge, old_vert_noKey);

      I8 e0_flip = -1;
      I8 e1_flip = -1;
      I8 e2_flip = -1;
      LO v1 = v1_old_face;
      LO v2 = v2_old_face;
      auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
      auto e0v1 = old_ev2v[old_face_e0*2 + 1];
      auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
      auto e1v1 = old_ev2v[old_face_e1*2 + 1];
      auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
      auto e2v1 = old_ev2v[old_face_e2*2 + 1];
      if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
        e0_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
      }
      if ((e1v0_old_face == v2) && (e1v1 == v1)) {
        e1_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
      }
      if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
        e2_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
      }

      Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
      Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
      //Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
      Real cx30 = old_vertCtrlPts[v1*dim + 0];
      Real cy30 = old_vertCtrlPts[v1*dim + 1];
      //Real cz30 = old_vertCtrlPts[v1*dim + 2];
      Real cx03 = old_vertCtrlPts[v2*dim + 0];
      Real cy03 = old_vertCtrlPts[v2*dim + 1];
      //Real cz03 = old_vertCtrlPts[v2*dim + 2];

      auto pts_per_edge = n_edge_pts;
      Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
      Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
      //Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
      Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
      Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        //auto tempz = cz10;
        cx10 = cx20;
        cy10 = cy20;
        //cz10 = cz20;
        cx20 = tempx;
        cy20 = tempy;
        //cz20 = tempz;
      }

      Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
      Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
      //Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
      Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        //auto tempz = cz21;
        cx21 = cx12;
        cy21 = cy12;
        //cz21 = cz12;
        cx12 = tempx;
        cy12 = tempy;
        //cz12 = tempz;
      }

      Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
      Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
      //Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
      Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        //auto tempz = cz02;
        cx02 = cx01;
        cy02 = cy01;
        //cz02 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        //cz01 = tempz;
      }

      Real cx11 = old_faceCtrlPts[old_face*dim + 0];
      Real cy11 = old_faceCtrlPts[old_face*dim + 1];
      //Real cz11 = old_faceCtrlPts[old_face*dim + 2];

      printf("ok1, new faces f0 f1 %d %d oldface %d\n", new_f0, new_f1, old_face);
      auto nodePts = cubic_face_xi_values
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
        auto p11_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
          cx10*B10_cube(nodePts[0], nodePts[1]) +
          cx20*B20_cube(nodePts[0], nodePts[1]) +
          cx30*B30_cube(nodePts[0], nodePts[1]) +
          cx21*B21_cube(nodePts[0], nodePts[1]) +
          cx12*B12_cube(nodePts[0], nodePts[1]) +
          cx03*B03_cube(nodePts[0], nodePts[1]) +
          cx02*B02_cube(nodePts[0], nodePts[1]) +
          cx01*B01_cube(nodePts[0], nodePts[1]) +
          cx11*B11_cube(nodePts[0], nodePts[1]);
        auto p11_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
          cy10*B10_cube(nodePts[0], nodePts[1]) +
          cy20*B20_cube(nodePts[0], nodePts[1]) +
          cy30*B30_cube(nodePts[0], nodePts[1]) +
          cy21*B21_cube(nodePts[0], nodePts[1]) +
          cy12*B12_cube(nodePts[0], nodePts[1]) +
          cy03*B03_cube(nodePts[0], nodePts[1]) +
          cy02*B02_cube(nodePts[0], nodePts[1]) +
          cy01*B01_cube(nodePts[0], nodePts[1]) +
          cy11*B11_cube(nodePts[0], nodePts[1]);

        //use these as interp pts to find ctrl pt in new face
        //inquire known vert and edge ctrl pts
        LO newface = new_f0;
        I8 newface_e0_flip = -1;
        I8 newface_e1_flip = -1;
        I8 newface_e2_flip = -1;
        LO newface_v0 = new_fv2v[newface*3 + 0];
        LO newface_v1 = new_fv2v[newface*3 + 1];
        LO newface_v2 = new_fv2v[newface*3 + 2];
        LO newface_e0 = new_fe2e[newface*3 + 0];
        LO newface_e1 = new_fe2e[newface*3 + 1];
        LO newface_e2 = new_fe2e[newface*3 + 2];
        auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
        auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
        auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
        auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
        auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
        auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
        if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
          newface_e0_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
        }
        if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
          newface_e1_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
        }
        if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
          newface_e2_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
        }

        printf("ok1\n");
        Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
        Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
        //Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
        Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
        Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
        //Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
        Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
        Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
        //Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];
        printf("ok2\n");

        Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
        Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
        //Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
        Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
        Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e0_flip > 0) {
          swap2(newface_cx10, newface_cx20);
          swap2(newface_cy10, newface_cy20);
          //swap2(newface_cz10, newface_cz20);
        }

        printf("ok3\n");
        Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
        Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
        //Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
        Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e1_flip > 0) {
          swap2(newface_cx21, newface_cx12);
          swap2(newface_cy21, newface_cy12);
          //swap2(newface_cz21, newface_cz12);
        }
        printf("ok4\n");

        Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
        Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
        //Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
        Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e2_flip > 0) {
          swap2(newface_cx01, newface_cx02);
          swap2(newface_cy01, newface_cy02);
          //swap2(newface_cz01, newface_cz02);
        }
        printf("ok5\n");

        auto xi_11 = xi_11_cube();
        Real newface_cx11 = (p11_x - newface_cx00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cx10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cx20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cx30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cx21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cx12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cx03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cx02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cx01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        Real newface_cy11 = (p11_y - newface_cy00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cy10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cy20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cy30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cy21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cy12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cy03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cy02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cy01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        if (old_face == 5)
        printf("for f0 newface %d p11 is (%f, %f) c11 is (%f,%f)\n",
               newface, p11_x, p11_y, newface_cx11, newface_cy11);
        face_ctrlPts[newface*n_face_pts*dim + 0] = newface_cx11;
        face_ctrlPts[newface*n_face_pts*dim + 1] = newface_cy11;
      }

      //for f1
      {
        //get the interp point
        auto p11_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
          cx10*B10_cube(nodePts[2], nodePts[3]) +
          cx20*B20_cube(nodePts[2], nodePts[3]) +
          cx30*B30_cube(nodePts[2], nodePts[3]) +
          cx21*B21_cube(nodePts[2], nodePts[3]) +
          cx12*B12_cube(nodePts[2], nodePts[3]) +
          cx03*B03_cube(nodePts[2], nodePts[3]) +
          cx02*B02_cube(nodePts[2], nodePts[3]) +
          cx01*B01_cube(nodePts[2], nodePts[3]) +
          cx11*B11_cube(nodePts[2], nodePts[3]);
        auto p11_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
          cy10*B10_cube(nodePts[2], nodePts[3]) +
          cy20*B20_cube(nodePts[2], nodePts[3]) +
          cy30*B30_cube(nodePts[2], nodePts[3]) +
          cy21*B21_cube(nodePts[2], nodePts[3]) +
          cy12*B12_cube(nodePts[2], nodePts[3]) +
          cy03*B03_cube(nodePts[2], nodePts[3]) +
          cy02*B02_cube(nodePts[2], nodePts[3]) +
          cy01*B01_cube(nodePts[2], nodePts[3]) +
          cy11*B11_cube(nodePts[2], nodePts[3]);

        //use these as interp pts to find ctrl pt in new face
        //inquire known vert and edge ctrl pts
        LO newface = new_f1;
        I8 newface_e0_flip = -1;
        I8 newface_e1_flip = -1;
        I8 newface_e2_flip = -1;
        LO newface_v0 = new_fv2v[newface*3 + 0];
        LO newface_v1 = new_fv2v[newface*3 + 1];
        LO newface_v2 = new_fv2v[newface*3 + 2];
        LO newface_e0 = new_fe2e[newface*3 + 0];
        LO newface_e1 = new_fe2e[newface*3 + 1];
        LO newface_e2 = new_fe2e[newface*3 + 2];
        auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
        auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
        auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
        auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
        auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
        auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
        if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
          newface_e0_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
        }
        if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
          newface_e1_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
        }
        if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
          newface_e2_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
        }

        Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
        Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
        //Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
        Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
        Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
        //Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
        Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
        Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
        //Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];

        Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
        Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
        //Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
        Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
        Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e0_flip > 0) {
          swap2(newface_cx10, newface_cx20);
          swap2(newface_cy10, newface_cy20);
          //swap2(newface_cz10, newface_cz20);
        }

        Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
        Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
        //Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
        Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e1_flip > 0) {
          swap2(newface_cx21, newface_cx12);
          swap2(newface_cy21, newface_cy12);
          //swap2(newface_cz21, newface_cz12);
        }

        Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
        Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
        //Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
        Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e2_flip > 0) {
          swap2(newface_cx01, newface_cx02);
          swap2(newface_cy01, newface_cy02);
          //swap2(newface_cz01, newface_cz02);
        }

        auto xi_11 = xi_11_cube();
        Real newface_cx11 = (p11_x - newface_cx00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cx10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cx20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cx30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cx21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cx12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cx03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cx02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cx01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        Real newface_cy11 = (p11_y - newface_cy00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cy10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cy20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cy30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cy21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cy12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cy03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cy02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cy01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        if (old_face == 5)
        printf("for f1 newface %d p11 is (%f, %f) c11 is (%f,%f)\n",
               newface, p11_x, p11_y, newface_cx11, newface_cy11);
        face_ctrlPts[newface*n_face_pts*dim + 0] = newface_cx11;
        face_ctrlPts[newface*n_face_pts*dim + 1] = newface_cy11;
      }
    }
    if ((end-start) == 3) {
      auto new_f0 = prods2new[start + 2];
      auto new_f1 = prods2new[start + 3];

      auto old_face = keys2old_faces[2*key + 1];
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
      printf("for old face %d , oldKeyEdge %d, found old no-key vert %d\n",
          old_face, old_key_edge, old_vert_noKey);

      I8 e0_flip = -1;
      I8 e1_flip = -1;
      I8 e2_flip = -1;
      LO v1 = v1_old_face;
      LO v2 = v2_old_face;
      auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
      auto e0v1 = old_ev2v[old_face_e0*2 + 1];
      auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
      auto e1v1 = old_ev2v[old_face_e1*2 + 1];
      auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
      auto e2v1 = old_ev2v[old_face_e2*2 + 1];
      if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
        e0_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
      }
      if ((e1v0_old_face == v2) && (e1v1 == v1)) {
        e1_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
      }
      if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
        e2_flip = 1;
      }
      else {
        OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
      }

      Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
      Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
      //Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
      Real cx30 = old_vertCtrlPts[v1*dim + 0];
      Real cy30 = old_vertCtrlPts[v1*dim + 1];
      //Real cz30 = old_vertCtrlPts[v1*dim + 2];
      Real cx03 = old_vertCtrlPts[v2*dim + 0];
      Real cy03 = old_vertCtrlPts[v2*dim + 1];
      //Real cz03 = old_vertCtrlPts[v2*dim + 2];

      auto pts_per_edge = n_edge_pts;
      Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
      Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
      //Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
      Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
      Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e0_flip > 0) {
        auto tempx = cx10;
        auto tempy = cy10;
        //auto tempz = cz10;
        cx10 = cx20;
        cy10 = cy20;
        //cz10 = cz20;
        cx20 = tempx;
        cy20 = tempy;
        //cz20 = tempz;
      }

      Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
      Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
      //Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
      Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e1_flip > 0) {
        auto tempx = cx21;
        auto tempy = cy21;
        //auto tempz = cz21;
        cx21 = cx12;
        cy21 = cy12;
        //cz21 = cz12;
        cx12 = tempx;
        cy12 = tempy;
        //cz12 = tempz;
      }

      Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
      Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
      //Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
      Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
      Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
      //Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
      if (e2_flip > 0) {
        auto tempx = cx02;
        auto tempy = cy02;
        //auto tempz = cz02;
        cx02 = cx01;
        cy02 = cy01;
        //cz02 = cz01;
        cx01 = tempx;
        cy01 = tempy;
        //cz01 = tempz;
      }

      Real cx11 = old_faceCtrlPts[old_face*dim + 0];
      Real cy11 = old_faceCtrlPts[old_face*dim + 1];
      //Real cz11 = old_faceCtrlPts[old_face*dim + 2];

      printf("ok3, new faces f2 f3 %d %d oldface %d\n", new_f0, new_f1, old_face);
      auto nodePts = cubic_face_xi_values
        (old_vert_noKey, v0_old_face, v1_old_face, v2_old_face, old_key_edge,
         old_face_e0, old_face_e1, old_face_e2, new_fv2v[new_f0*3 + 0], 
         new_fv2v[new_f0*3 + 1], new_fv2v[new_f0*3 + 2], new_fv2v[new_f1*3 + 0], 
         new_fv2v[new_f1*3 + 1], new_fv2v[new_f1*3 + 2],
         old_verts2new_verts[v0_old_face],
         old_verts2new_verts[v1_old_face], old_verts2new_verts[v2_old_face]);
      printf("ok4\n");

      //for f2
      {
        //get the interp point
        auto p11_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
          cx10*B10_cube(nodePts[0], nodePts[1]) +
          cx20*B20_cube(nodePts[0], nodePts[1]) +
          cx30*B30_cube(nodePts[0], nodePts[1]) +
          cx21*B21_cube(nodePts[0], nodePts[1]) +
          cx12*B12_cube(nodePts[0], nodePts[1]) +
          cx03*B03_cube(nodePts[0], nodePts[1]) +
          cx02*B02_cube(nodePts[0], nodePts[1]) +
          cx01*B01_cube(nodePts[0], nodePts[1]) +
          cx11*B11_cube(nodePts[0], nodePts[1]);
        auto p11_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
          cy10*B10_cube(nodePts[0], nodePts[1]) +
          cy20*B20_cube(nodePts[0], nodePts[1]) +
          cy30*B30_cube(nodePts[0], nodePts[1]) +
          cy21*B21_cube(nodePts[0], nodePts[1]) +
          cy12*B12_cube(nodePts[0], nodePts[1]) +
          cy03*B03_cube(nodePts[0], nodePts[1]) +
          cy02*B02_cube(nodePts[0], nodePts[1]) +
          cy01*B01_cube(nodePts[0], nodePts[1]) +
          cy11*B11_cube(nodePts[0], nodePts[1]);

        //use these as interp pts to find ctrl pt in new face
        //inquire known vert and edge ctrl pts
        LO newface = new_f0;
        LO newface_e0 = new_fe2e[newface*3 + 0];
        LO newface_e1 = new_fe2e[newface*3 + 1];
        LO newface_e2 = new_fe2e[newface*3 + 2];
        I8 newface_e0_flip = -1;
        I8 newface_e1_flip = -1;
        I8 newface_e2_flip = -1;
        LO newface_v0 = new_fv2v[newface*3 + 0];
        LO newface_v1 = new_fv2v[newface*3 + 1];
        LO newface_v2 = new_fv2v[newface*3 + 2];
        auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
        auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
        auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
        auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
        auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
        auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
        if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
          newface_e0_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
        }
        if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
          newface_e1_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
        }
        if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
          newface_e2_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
        }

        printf("ok1\n");
        Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
        Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
        //Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
        Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
        Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
        //Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
        Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
        Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
        //Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];
        printf("ok2\n");

        Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
        Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
        //Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
        Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
        Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e0_flip > 0) {
          swap2(newface_cx10, newface_cx20);
          swap2(newface_cy10, newface_cy20);
          //swap2(newface_cz10, newface_cz20);
        }

        printf("ok3\n");
        Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
        Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
        //Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
        Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e1_flip > 0) {
          swap2(newface_cx21, newface_cx12);
          swap2(newface_cy21, newface_cy12);
          //swap2(newface_cz21, newface_cz12);
        }
        printf("ok4\n");

        Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
        Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
        //Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
        Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e2_flip > 0) {
          swap2(newface_cx01, newface_cx02);
          swap2(newface_cy01, newface_cy02);
          //swap2(newface_cz01, newface_cz02);
        }
        printf("ok5\n");

        auto xi_11 = xi_11_cube();
        Real newface_cx11 = (p11_x - newface_cx00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cx10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cx20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cx30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cx21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cx12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cx03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cx02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cx01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        Real newface_cy11 = (p11_y - newface_cy00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cy10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cy20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cy30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cy21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cy12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cy03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cy02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cy01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        if (old_face == 5)
          printf("for f2 newface %d p11 is (%f, %f) c11 is (%f,%f)\n",
                 newface, p11_x, p11_y, newface_cx11, newface_cy11);
        face_ctrlPts[newface*n_face_pts*dim + 0] = newface_cx11;
        face_ctrlPts[newface*n_face_pts*dim + 1] = newface_cy11;
      }

      //for f3
      {
        //get the interp point
        auto p11_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
          cx10*B10_cube(nodePts[2], nodePts[3]) +
          cx20*B20_cube(nodePts[2], nodePts[3]) +
          cx30*B30_cube(nodePts[2], nodePts[3]) +
          cx21*B21_cube(nodePts[2], nodePts[3]) +
          cx12*B12_cube(nodePts[2], nodePts[3]) +
          cx03*B03_cube(nodePts[2], nodePts[3]) +
          cx02*B02_cube(nodePts[2], nodePts[3]) +
          cx01*B01_cube(nodePts[2], nodePts[3]) +
          cx11*B11_cube(nodePts[2], nodePts[3]);
        auto p11_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
          cy10*B10_cube(nodePts[2], nodePts[3]) +
          cy20*B20_cube(nodePts[2], nodePts[3]) +
          cy30*B30_cube(nodePts[2], nodePts[3]) +
          cy21*B21_cube(nodePts[2], nodePts[3]) +
          cy12*B12_cube(nodePts[2], nodePts[3]) +
          cy03*B03_cube(nodePts[2], nodePts[3]) +
          cy02*B02_cube(nodePts[2], nodePts[3]) +
          cy01*B01_cube(nodePts[2], nodePts[3]) +
          cy11*B11_cube(nodePts[2], nodePts[3]);

        //use these as interp pts to find ctrl pt in new face
        //inquire known vert and edge ctrl pts
        LO newface = new_f1;
        I8 newface_e0_flip = -1;
        I8 newface_e1_flip = -1;
        I8 newface_e2_flip = -1;
        LO newface_v0 = new_fv2v[newface*3 + 0];
        LO newface_v1 = new_fv2v[newface*3 + 1];
        LO newface_v2 = new_fv2v[newface*3 + 2];
        LO newface_e0 = new_fe2e[newface*3 + 0];
        LO newface_e1 = new_fe2e[newface*3 + 1];
        LO newface_e2 = new_fe2e[newface*3 + 2];
        auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
        auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
        auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
        auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
        auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
        auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
        if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
          newface_e0_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
        }
        if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
          newface_e1_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
        }
        if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
          newface_e2_flip = 1;
        }
        else {
          OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
        }

        Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
        Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
        //Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
        Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
        Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
        //Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
        Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
        Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
        //Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];

        Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
        Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
        //Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
        Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
        Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e0_flip > 0) {
          swap2(newface_cx10, newface_cx20);
          swap2(newface_cy10, newface_cy20);
          //swap2(newface_cz10, newface_cz20);
        }

        Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
        Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
        //Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
        Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e1_flip > 0) {
          swap2(newface_cx21, newface_cx12);
          swap2(newface_cy21, newface_cy12);
          //swap2(newface_cz21, newface_cz12);
        }

        Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
        Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
        //Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
        Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
        Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
        //Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
        if (newface_e2_flip > 0) {
          swap2(newface_cx01, newface_cx02);
          swap2(newface_cy01, newface_cy02);
          //swap2(newface_cz01, newface_cz02);
        }

        auto xi_11 = xi_11_cube();
        Real newface_cx11 = (p11_x - newface_cx00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cx10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cx20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cx30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cx21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cx12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cx03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cx02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cx01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        Real newface_cy11 = (p11_y - newface_cy00*B00_cube(xi_11[0], xi_11[1]) -
            newface_cy10*B10_cube(xi_11[0], xi_11[1]) -
            newface_cy20*B20_cube(xi_11[0], xi_11[1]) -
            newface_cy30*B30_cube(xi_11[0], xi_11[1]) -
            newface_cy21*B21_cube(xi_11[0], xi_11[1]) -
            newface_cy12*B12_cube(xi_11[0], xi_11[1]) -
            newface_cy03*B03_cube(xi_11[0], xi_11[1]) -
            newface_cy02*B02_cube(xi_11[0], xi_11[1]) -
            newface_cy01*B01_cube(xi_11[0], xi_11[1]))/B11_cube(xi_11[0], xi_11[1]);
        if (old_face == 5)
        printf("for f3 newface %d p11 is (%f, %f) c11 is (%f,%f)\n",
               newface, p11_x, p11_y, newface_cx11, newface_cy11);
        face_ctrlPts[newface*n_face_pts*dim + 0] = newface_cx11;
        face_ctrlPts[newface*n_face_pts*dim + 1] = newface_cy11;
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
        OMEGA_H_CHECK(face_ctrlPts[new_face*n_face_pts*dim + d] = INT8_MAX);
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
