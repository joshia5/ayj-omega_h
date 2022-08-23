#include "Omega_h_for.hpp"
#include "Omega_h_curve_validity_3d.hpp"
#include "Omega_h_coarsen.hpp"
#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_coarsen_invalidities.hpp"

#include <iostream>

namespace Omega_h {

LOs coarsen_invalidities_2d(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  LO const mesh_dim = 2;
  OMEGA_H_CHECK(mesh->dim() == 2);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_elem_verts();
  auto v2c = mesh->ask_up(VERT, mesh_dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto invalidities = Write<LO>(ncands * 2, -1);
  auto coords = mesh->coords();
  auto e2e = edges_across_tris(mesh->ask_down(2, 1), mesh->ask_up(1, 2));
  auto e2ee = e2e.a2ab;
  auto ee2e = e2e.ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      auto eev_onto = 1 - eev_col;
      auto v_onto = ev2v[e * 2 + eev_onto];
      LO max_invalid = -1;
      for (auto vc = v2vc[v_col]; vc < v2vc[v_col + 1]; ++vc) {
        auto c = vc2c[vc];
        auto vc_code = vc_codes[vc];
        auto ccv_col = code_which_down(vc_code);
        auto ccv2v = gather_verts<2 + 1>(cv2v, c);
        bool will_die = false;
        for (auto ccv = 0; ccv < (2 + 1); ++ccv) {
          if ((ccv != ccv_col) && (ccv2v[ccv] == v_onto)) {
            will_die = true;
            break;
          }
        }
        if (will_die) continue;
        OMEGA_H_CHECK(0 <= ccv_col && ccv_col < mesh_dim + 1);
        ccv2v[ccv_col] = v_onto;  // vertices of new cell

        Few<LO, 2> same_edges;
        for (LO init = 0; init < 2; ++init) {
          same_edges[init] = -1; //can be max 2
        }
        Few<LO, 2+1> is_newTri_edge_flip;
        for (LO init = 0; init < 3; ++init) {
          is_newTri_edge_flip[init] = -2;
        }
        Few<LO, 2+1> newTri_edge;
        for (LO init = 0; init < 3; ++init) {
          newTri_edge[init] = -1;
        }
        LO count_same_edge = 0;
        auto v0_f = ccv2v[0];
        auto v1_f = ccv2v[1];
        auto v2_f = ccv2v[2];
        for (auto ee = e2ee[e]; ee < e2ee[e + 1]; ++ee) {
          if (count_same_edge >= 2) break;
          auto adj_e = ee2e[ee];
          //printf("e %d adj_e %d\n", e , adj_e);
          auto v0 = ev2v[adj_e*2];
          auto v1 = ev2v[adj_e*2 + 1];

          //check first edge
          if (((v0 == v0_f) && (v1 == v1_f)) ||
              ((v0 == v1_f) && (v1 == v0_f))) {
            if ((count_same_edge == 0) || 
                ((count_same_edge == 1) && 
                 (same_edges[count_same_edge-1] != adj_e))) {
              same_edges[count_same_edge] = adj_e;
              newTri_edge[0] = adj_e;
              if ((v0 == v1_f) && (v1 == v0_f)) {
                is_newTri_edge_flip[0] = 1;
              }
              else {
                is_newTri_edge_flip[0] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 2nd edge
          else if (((v0 == v1_f) && (v1 == v2_f)) ||
              ((v0 == v2_f) && (v1 == v1_f))) {
            if ((count_same_edge == 0) || 
                ((count_same_edge == 1) && 
                 (same_edges[count_same_edge-1] != adj_e))) {
              same_edges[count_same_edge] = adj_e;
              newTri_edge[1] = adj_e;
              if ((v0 == v2_f) && (v1 == v1_f)) {
                is_newTri_edge_flip[1] = 1;
              }
              else {
                is_newTri_edge_flip[1] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 3nd edge
          else if (((v0 == v2_f) && (v1 == v0_f)) ||
              ((v0 == v0_f) && (v1 == v2_f))) {
            if ((count_same_edge == 0) || 
                ((count_same_edge == 1) && 
                 (same_edges[count_same_edge-1] != adj_e))) {
              same_edges[count_same_edge] = adj_e;
              newTri_edge[2] = adj_e;
              if ((v0 == v0_f) && (v1 == v2_f)) {
                is_newTri_edge_flip[2] = 1;
              }
              else {
                is_newTri_edge_flip[2] = -1;
              }
              ++count_same_edge;
            }
          }
          else {
            for (auto ee2 = e2ee[adj_e]; ee2 < e2ee[adj_e + 1]; ++ee2) {
              if (count_same_edge >= 2) break;
              auto adj_e2 = ee2e[ee2];
              auto v0_e2 = ev2v[adj_e2*2];
              auto v1_e2 = ev2v[adj_e2*2 + 1];

              //check first edge
              if (((v0_e2  == v0_f) && (v1_e2  == v1_f)) ||
                  ((v0_e2  == v1_f) && (v1_e2  == v0_f))) {
                if ((count_same_edge == 0) || 
                    ((count_same_edge == 1) && 
                     (same_edges[count_same_edge-1] != adj_e2))) {
                  same_edges[count_same_edge] = adj_e2;
                  newTri_edge[0] = adj_e2;
                  if ((v0_e2  == v1_f) && (v1_e2  == v0_f)) {
                    is_newTri_edge_flip[0] = 1;
                  }
                  else {
                    is_newTri_edge_flip[0] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 2nd edge
              else if (((v0_e2 == v1_f) && (v1_e2 == v2_f)) ||
                  ((v0_e2 == v2_f) && (v1_e2 == v1_f))) {
                if ((count_same_edge == 0) || 
                    ((count_same_edge == 1) && 
                     (same_edges[count_same_edge-1] != adj_e2))) {
                  same_edges[count_same_edge] = adj_e2;
                  newTri_edge[1] = adj_e2;
                  if ((v0_e2 == v2_f) && (v1_e2 == v1_f)) {
                    is_newTri_edge_flip[1] = 1;
                  }
                  else {
                    is_newTri_edge_flip[1] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 3nd edge
              else if (((v0_e2 == v2_f) && (v1_e2 == v0_f)) ||
                  ((v0_e2 == v0_f) && (v1_e2 == v2_f))) {
                if ((count_same_edge == 0) || 
                    ((count_same_edge == 1) && 
                     (same_edges[count_same_edge-1] != adj_e2))) {
                  same_edges[count_same_edge] = adj_e2;
                  newTri_edge[2] = adj_e2;
                  if ((v0_e2 == v0_f) && (v1_e2 == v2_f)) {
                    is_newTri_edge_flip[2] = 1;
                  }
                  else {
                    is_newTri_edge_flip[2] = -1;
                  }
                  ++count_same_edge;
                }
              }
              else {
                //printf("adje2 %d not same edge \n", adj_e2);
              }
            }
          }

        }
        for (LO ee = 0; ee < 2; ++ee) {
          //printf("cand %d, e %d, found same edge %d with newTri %d %d %d vcol %d vOnto %d\n",
          //cand, e, same_edges[ee], v0_f, v1_f, v2_f, v_col, v_onto);
        }
        for (LO ee = 0; ee < 3; ++ee) {
          //printf("newtri edge %d is %d\n", ee, newTri_edge[ee]);
        }
        Few<Real, 10*2> tri_pts;//ntri_pts*dim=20
        for (LO j = 0; j < 3; ++j) {
          auto p = get_vector<2>(vertCtrlPts, ccv2v[j]);
          for (LO k = 0; k < mesh_dim; ++k) {
            tri_pts[j*mesh_dim + k] = p[k];
          }
        }
        Few<Real, 3*2> thirdOfEdge;
        for (I8 d = 0; d < mesh_dim; ++d) {
          thirdOfEdge[0*mesh_dim + d] = 
            1.0/3.0*(tri_pts[1*mesh_dim + d] - tri_pts[0*mesh_dim + d]);
          thirdOfEdge[1*mesh_dim + d] = 
            1.0/3.0*(tri_pts[2*mesh_dim + d] - tri_pts[1*mesh_dim + d]);
          thirdOfEdge[2*mesh_dim + d] = 
            1.0/3.0*(tri_pts[0*mesh_dim + d] - tri_pts[2*mesh_dim + d]);
        }
        for (LO j = 0; j < 3; ++j) {
          LO index = 3;
          if (is_newTri_edge_flip[j] == -1) {
            for (I8 d = 0; d < mesh_dim; ++d) {
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
                edgeCtrlPts[newTri_edge[j]*n_edge_pts*mesh_dim + d];
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                edgeCtrlPts[newTri_edge[j]*n_edge_pts*mesh_dim + mesh_dim + d];
            }
          }
          else if (is_newTri_edge_flip[j] == 1) {
            for (I8 d = 0; d < mesh_dim; ++d) {
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
                edgeCtrlPts[newTri_edge[j]*n_edge_pts*mesh_dim + mesh_dim + d];
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                edgeCtrlPts[newTri_edge[j]*n_edge_pts*mesh_dim + d];
            }
          }
          else {
            for (I8 d = 0; d < mesh_dim; ++d) {
              assert (is_newTri_edge_flip[j] == -2);
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] = 
                tri_pts[j*mesh_dim + d] + thirdOfEdge[j*mesh_dim + d];
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                tri_pts[j*mesh_dim + d] + 2*thirdOfEdge[j*mesh_dim + d];
            }
          }
        }
        //query the face's ctrl pt and store
        //TODO triPt using blending?
        for (I8 d = 0; d < mesh_dim; ++d) {
          LO index = 9;
          tri_pts[index*mesh_dim + d] = 1.0/3.0*(
              tri_pts[0*mesh_dim + d] + tri_pts[1*mesh_dim + d] + 
              tri_pts[2*mesh_dim + d]);
        }
        auto nodes_det = getTriJacDetNodes<15, mesh_dim>(3, tri_pts);
        auto is_invalid = checkMinJacDet<15>(nodes_det, 3);
        //printf("cand %d eev_col %d, vc = %d is invalid %d\n", cand, eev_col, vc, is_invalid);

        max_invalid = max2(max_invalid, is_invalid);
      }
      invalidities[cand * 2 + eev_col] = max_invalid;
    }
  };
  parallel_for(ncands, f, "coarsen_invalidities");
  auto out_invalid = LOs(invalidities);
  return mesh->sync_subset_array(EDGE, out_invalid, cands2edges, -1, 2);
}

LOs coarsen_invalidities_3d(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  LO const mesh_dim = 3;
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_elem_verts();
  auto v2c = mesh->ask_up(VERT, mesh_dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto invalidities = Write<LO>(ncands * 2, -1);
  auto coords = mesh->coords();
  auto e2e = edges_across_tets(mesh->ask_down(3, 1), mesh->ask_up(1, 3));
  auto e2ee = e2e.a2ab;
  auto ee2e = e2e.ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      auto eev_onto = 1 - eev_col;
      auto v_onto = ev2v[e * 2 + eev_onto];
      LO max_invalid = -1;
      for (auto vc = v2vc[v_col]; vc < v2vc[v_col + 1]; ++vc) {
        auto c = vc2c[vc];
        auto vc_code = vc_codes[vc];
        auto ccv_col = code_which_down(vc_code);
        auto ccv2v = gather_verts<3 + 1>(cv2v, c);
        bool will_die = false;
        for (auto ccv = 0; ccv < (mesh_dim + 1); ++ccv) {
          if ((ccv != ccv_col) && (ccv2v[ccv] == v_onto)) {
            will_die = true;
            break;
          }
        }
        if (will_die) continue;
        OMEGA_H_CHECK(0 <= ccv_col && ccv_col < mesh_dim + 1);
        ccv2v[ccv_col] = v_onto;  // vertices of new cell

        Few<LO, 5> same_edges; //can be max 5
        for (LO init = 0; init < 5; ++init) {
          same_edges[init] = -1;
        }
        Few<LO, 6> is_newTet_edge_flip;
        for (LO init = 0; init < 6; ++init) {
          is_newTet_edge_flip[init] = -2;
        }
        Few<LO, 6> newTet_edge;
        for (LO init = 0; init < 6; ++init) {
          newTet_edge[init] = -1;
        }
        LO count_same_edge = 0;
        auto v0_r = ccv2v[0];
        auto v1_r = ccv2v[1];
        auto v2_r = ccv2v[2];
        auto v3_r = ccv2v[3];
        for (auto ee = e2ee[e]; ee < e2ee[e + 1]; ++ee) {
          if (count_same_edge >= 5) break;
          auto adj_e = ee2e[ee];
          //printf("e %d adj_e %d\n", e , adj_e);
          auto v0 = ev2v[adj_e*2];
          auto v1 = ev2v[adj_e*2 + 1];

          //check first edge
          if (((v0 == v0_r) && (v1 == v1_r)) ||
              ((v0 == v1_r) && (v1 == v0_r))) {
            if (newTet_edge[0] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[0] = adj_e;
              if ((v0 == v1_r) && (v1 == v0_r)) {
                is_newTet_edge_flip[0] = 1;
              }
              else {
                is_newTet_edge_flip[0] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 2nd edge
          else if (((v0 == v1_r) && (v1 == v2_r)) ||
              ((v0 == v2_r) && (v1 == v1_r))) {
            if (newTet_edge[1] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[1] = adj_e;
              if ((v0 == v2_r) && (v1 == v1_r)) {
                is_newTet_edge_flip[1] = 1;
              }
              else {
                is_newTet_edge_flip[1] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 3nd edge
          else if (((v0 == v2_r) && (v1 == v0_r)) ||
              ((v0 == v0_r) && (v1 == v2_r))) {
            if (newTet_edge[2] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[2] = adj_e;
              if ((v0 == v0_r) && (v1 == v2_r)) {
                is_newTet_edge_flip[2] = 1;
              }
              else {
                is_newTet_edge_flip[2] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 4th edge
          else if (((v0 == v0_r) && (v1 == v3_r)) ||
              ((v0 == v3_r) && (v1 == v0_r))) {
            if (newTet_edge[3] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[3] = adj_e;
              if ((v0 == v3_r) && (v1 == v0_r)) {
                is_newTet_edge_flip[3] = 1;
              }
              else {
                is_newTet_edge_flip[3] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 5th edge
          else if (((v0 == v1_r) && (v1 == v3_r)) ||
              ((v0 == v3_r) && (v1 == v1_r))) {
            if (newTet_edge[4] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[4] = adj_e;
              if ((v0 == v3_r) && (v1 == v1_r)) {
                is_newTet_edge_flip[4] = 1;
              }
              else {
                is_newTet_edge_flip[4] = -1;
              }
              ++count_same_edge;
            }
          }
          //check 6th edge
          else if (((v0 == v2_r) && (v1 == v3_r)) ||
              ((v0 == v3_r) && (v1 == v2_r))) {
            if (newTet_edge[5] == -1) {
              same_edges[count_same_edge] = adj_e;
              newTet_edge[5] = adj_e;
              if ((v0 == v3_r) && (v1 == v2_r)) {
                is_newTet_edge_flip[5] = 1;
              }
              else {
                is_newTet_edge_flip[5] = -1;
              }
              ++count_same_edge;
            }
          }
          else {
            for (auto ee2 = e2ee[adj_e]; ee2 < e2ee[adj_e + 1]; ++ee2) {
              if (count_same_edge >= 5) break;
              auto adj_e2 = ee2e[ee2];
              auto v0_e2 = ev2v[adj_e2*2];
              auto v1_e2 = ev2v[adj_e2*2 + 1];

              //check first edge
              if (((v0_e2 == v0_r) && (v1_e2 == v1_r)) ||
                  ((v0_e2 == v1_r) && (v1_e2 == v0_r))) {
                if (newTet_edge[0] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[0] = adj_e2;
                  if ((v0_e2 == v1_r) && (v1_e2 == v0_r)) {
                    is_newTet_edge_flip[0] = 1;
                  }
                  else {
                    is_newTet_edge_flip[0] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 2nd edge
              else if (((v0_e2 == v1_r) && (v1_e2 == v2_r)) ||
                  ((v0_e2 == v2_r) && (v1_e2 == v1_r))) {
                if (newTet_edge[1] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[1] = adj_e2;
                  if ((v0_e2 == v2_r) && (v1_e2 == v1_r)) {
                    is_newTet_edge_flip[1] = 1;
                  }
                  else {
                    is_newTet_edge_flip[1] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 3nd edge
              else if (((v0_e2 == v2_r) && (v1_e2 == v0_r)) ||
                  ((v0_e2 == v0_r) && (v1_e2 == v2_r))) {
                if (newTet_edge[2] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[2] = adj_e2;
                  if ((v0_e2 == v0_r) && (v1_e2 == v2_r)) {
                    is_newTet_edge_flip[2] = 1;
                  }
                  else {
                    is_newTet_edge_flip[2] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 4th edge
              else if (((v0_e2 == v0_r) && (v1_e2 == v3_r)) ||
                  ((v0_e2 == v3_r) && (v1_e2 == v0_r))) {
                if (newTet_edge[3] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[3] = adj_e2;
                  if ((v0_e2 == v3_r) && (v1_e2 == v0_r)) {
                    is_newTet_edge_flip[3] = 1;
                  }
                  else {
                    is_newTet_edge_flip[3] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 5th edge
              else if (((v0_e2 == v1_r) && (v1_e2 == v3_r)) ||
                  ((v0_e2 == v3_r) && (v1_e2 == v1_r))) {
                if (newTet_edge[4] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[4] = adj_e2;
                  if ((v0_e2 == v3_r) && (v1_e2 == v1_r)) {
                    is_newTet_edge_flip[4] = 1;
                  }
                  else {
                    is_newTet_edge_flip[4] = -1;
                  }
                  ++count_same_edge;
                }
              }
              //check 6th edge
              else if (((v0_e2 == v2_r) && (v1_e2 == v3_r)) ||
                  ((v0_e2 == v3_r) && (v1_e2 == v2_r))) {
                if (newTet_edge[5] == -1) {
                  same_edges[count_same_edge] = adj_e2;
                  newTet_edge[5] = adj_e2;
                  if ((v0_e2 == v3_r) && (v1_e2 == v2_r)) {
                    is_newTet_edge_flip[5] = 1;
                  }
                  else {
                    is_newTet_edge_flip[5] = -1;
                  }
                  ++count_same_edge;
                }
              }
              else {
                //printf("adje2 %d not same edge \n", adj_e2);
              }
            }
          }
        }
        for (LO ee = 0; ee < 2; ++ee) {
          //printf("cand %d, e %d, found same edge %d with newTri %d %d %d vcol %d vOnto %d\n",
          //cand, e, same_edges[ee], v0_f, v1_f, v2_f, v_col, v_onto);
        }
        for (LO ee = 0; ee < 3; ++ee) {
          //printf("newtri edge %d is %d\n", ee, newTri_edge[ee]);
        }
        Few<Real, 20*3> tet_pts;//ntri_pts*dim=20
        for (LO j = 0; j < 4; ++j) {
          auto p = get_vector<3>(vertCtrlPts, ccv2v[j]);
          for (LO k = 0; k < mesh_dim; ++k) {
            tet_pts[j*mesh_dim + k] = p[k];
          }
        }
        Few<Real, 6*3> thirdOfEdge;
        for (I8 d = 0; d < mesh_dim; ++d) {
          thirdOfEdge[0*mesh_dim + d] = 
            1.0/3.0*(tet_pts[1*mesh_dim + d] - tet_pts[0*mesh_dim + d]);
          thirdOfEdge[1*mesh_dim + d] = 
            1.0/3.0*(tet_pts[2*mesh_dim + d] - tet_pts[1*mesh_dim + d]);
          thirdOfEdge[2*mesh_dim + d] = 
            1.0/3.0*(tet_pts[0*mesh_dim + d] - tet_pts[2*mesh_dim + d]);
          thirdOfEdge[3*mesh_dim + d] = 
            1.0/3.0*(tet_pts[3*mesh_dim + d] - tet_pts[0*mesh_dim + d]);
          thirdOfEdge[4*mesh_dim + d] = 
            1.0/3.0*(tet_pts[3*mesh_dim + d] - tet_pts[1*mesh_dim + d]);
          thirdOfEdge[5*mesh_dim + d] = 
            1.0/3.0*(tet_pts[3*mesh_dim + d] - tet_pts[2*mesh_dim + d]);
        }
        for (LO j = 0; j < 6; ++j) {
          LO index = 4;
          if (is_newTet_edge_flip[j] == -1) {
            //same edge not flipped
            for (I8 d = 0; d < mesh_dim; ++d) {
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
                edgeCtrlPts[newTet_edge[j]*n_edge_pts*mesh_dim + d];
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                edgeCtrlPts[newTet_edge[j]*n_edge_pts*mesh_dim + mesh_dim + d];
            }
          }
          else if (is_newTet_edge_flip[j] == 1) {
            //same edge flipped
            for (I8 d = 0; d < mesh_dim; ++d) {
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
                edgeCtrlPts[newTet_edge[j]*n_edge_pts*mesh_dim + mesh_dim + d];
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                edgeCtrlPts[newTet_edge[j]*n_edge_pts*mesh_dim + d];
            }
          }
          else {
            //new edge (straight sided)
            for (I8 d = 0; d < mesh_dim; ++d) {
              assert (is_newTet_edge_flip[j] == -2);
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] = 
                tet_pts[(j%3)*mesh_dim + d] + 3*thirdOfEdge[j*mesh_dim + d]*xi_1_cube();
              //start from start of edge so j%2
              tet_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
                tet_pts[(j%3)*mesh_dim + d] + 3*thirdOfEdge[j*mesh_dim + d]*xi_2_cube();
            }
          }
        }
        //query the face's ctrl pt and store
        //TODO triPt using blending?
        //TODO find same faces using adjacency search
        for (I8 d = 0; d < mesh_dim; ++d) {
          LO index = 16;
          tet_pts[index*mesh_dim + d] = 1.0/3.0*(
              tet_pts[0*mesh_dim + d] + tet_pts[1*mesh_dim + d] + 
              tet_pts[2*mesh_dim + d]);
          tet_pts[(index+1)*mesh_dim + d] = 1.0/3.0*(
              tet_pts[0*mesh_dim + d] + tet_pts[1*mesh_dim + d] + 
              tet_pts[3*mesh_dim + d]);
          tet_pts[(index+2)*mesh_dim + d] = 1.0/3.0*(
              tet_pts[1*mesh_dim + d] + tet_pts[2*mesh_dim + d] + 
              tet_pts[3*mesh_dim + d]);
          tet_pts[(index+3)*mesh_dim + d] = 1.0/3.0*(
              tet_pts[0*mesh_dim + d] + tet_pts[2*mesh_dim + d] + 
              tet_pts[3*mesh_dim + d]);
        }
        auto nodes_det = getTetJacDetNodes<84>(3, tet_pts);
        auto is_invalid = checkMinJacDet_3d(nodes_det, 3);

        max_invalid = max2(max_invalid, is_invalid);
      }
      invalidities[cand * 2 + eev_col] = max_invalid;
    }
  };
  parallel_for(ncands, f, "coarsen_invalidities");
  auto out_invalid = LOs(invalidities);
  return mesh->sync_subset_array(EDGE, out_invalid, cands2edges, -1, 2);
}

LOs coarsen_invalidities_new_mesh(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes, Mesh* new_mesh,
    LOs const old_verts2new_verts, Read<I8> const verts_are_keys) {
  LO const mesh_dim = mesh->dim();
  OMEGA_H_CHECK(mesh_dim == 3);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto ncands = cands2edges.size();
  auto invalidities = Write<LO>(ncands*2, -1);

  auto v2vc_new = (new_mesh->ask_up(VERT, mesh_dim)).a2ab;
  auto vc2c_new = (new_mesh->ask_up(VERT, mesh_dim)).ab2b;

  auto const new_ev2v = new_mesh->ask_down(1, 0).ab2b;
  auto const new_rv2v = new_mesh->ask_down(3, 0).ab2b;
  auto const new_re2e = new_mesh->ask_down(3, 1).ab2b;
  auto const new_rf2f = new_mesh->ask_down(3, 2).ab2b;
  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const faceCtrlPts = new_mesh->get_ctrlPts(2);

  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e*2 + eev_col];
      //printf("vcol is key %d\n", verts_are_keys[v_col]);
      if (verts_are_keys[v_col] < 1) continue;
      auto eev_onto = 1 - eev_col;
      auto v_onto = ev2v[e*2 + eev_onto];
      auto v_onto_new = old_verts2new_verts[v_onto];
      if (v_onto_new == -1) continue;
      printf("edge %d v_onto_new %d\n", e, v_onto_new);
      LO max_invalid = -1;
      for (LO vc = v2vc_new[v_onto_new]; vc < v2vc_new[v_onto_new+1]; ++vc) {
        auto c = vc2c_new[vc];
        printf("new tet %d\n", c);
        Few<Real, 60> tet_pts = collect_tet_pts(3, c, new_ev2v, new_rv2v,
          vertCtrlPts, edgeCtrlPts, faceCtrlPts, new_re2e, new_rf2f);

        auto nodes_det = getTetJacDetNodes<84>(3, tet_pts);
        auto is_invalid = checkMinJacDet_3d(nodes_det, 3);

        max_invalid = max2(max_invalid, is_invalid);
      }
      invalidities[cand*2 + eev_col] = max_invalid;
    }
  };
  parallel_for(ncands, f, "coarsen_invalidities");
  auto out_invalid = LOs(invalidities);
  return mesh->sync_subset_array(EDGE, out_invalid, cands2edges, -1, 2);
}

LOs coarsen_invalidities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  if (mesh->dim() == 3) {
    return coarsen_invalidities_3d(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2) {
    return coarsen_invalidities_2d(mesh, cands2edges, cand_codes);
  }
  OMEGA_H_NORETURN(LOs());
}

Read<I8> filter_coarsen_invalids(
    Read<I8> cand_codes, LOs cand_invalids, LO is_invalid) {
  auto keep_dirs = each_eq_to(cand_invalids, is_invalid);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

}  // end namespace Omega_h
