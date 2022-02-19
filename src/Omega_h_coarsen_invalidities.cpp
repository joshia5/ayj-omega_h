#include "Omega_h_coarsen.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_metric.hpp"

#include <iostream>

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
LOs coarsen_invalidities_tmpl(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->dim() == mesh_dim);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_elem_verts();
  auto v2c = mesh->ask_up(VERT, mesh_dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto invalidities = Write<LO>(ncands * 2, -1);
  auto coords = mesh->coords();

  auto e2f = mesh->ask_up(1, 2);
  auto f2e = mesh->ask_down(2, 1);
  auto e2e = edges_across_tris(f2e, e2f);
  auto e2ee = e2e.a2ab;
  auto ee2e = e2e.ab2b;

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
        auto ccv2v = gather_verts<mesh_dim + 1>(cv2v, c);
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

        //TODO use these and interpolate edges to calculate ctrl pts and give
        //them to validity check
        //so there will be 2 functions
        //1. construct_elemNodes_from_adj
        //2. construct_elemNodes_from_verts
        //3. only the new edge is straight, others are same as old
        Few<LO, 2> same_edges = {-1, -1}; //can be max 2
        Few<LO, 3> is_newTri_edge_flip = {-2, -2, -2};
        Few<LO, 3> newTri_edge = {-1, -1, -1};
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

          //TODO this might also be a good place to check for flip
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
        auto vertCtrlPts = mesh->get_ctrlPts(0);
        auto edgeCtrlPts = mesh->get_ctrlPts(1);
        auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
        Few<Real, 20> tri_pts;//ntri_pts*dim=20
        for (LO j = 0; j < 3; ++j) {
          auto p = get_vector<2>(vertCtrlPts, ccv2v[j]);
          for (LO k = 0; k < mesh_dim; ++k) {
            tri_pts[j*mesh_dim + k] = p[k];
          }
        }
        Few<Real, 3*mesh_dim> thirdOfEdge;
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
        for (I8 d = 0; d < mesh_dim; ++d) {
          LO index = 9;
          tri_pts[index*mesh_dim + d] = 1.0/3.0*(
              tri_pts[0*mesh_dim + d] + tri_pts[1*mesh_dim + d] + 
              tri_pts[2*mesh_dim + d]);
        }
        auto nodes_det = getTriJacDetNodes<15>(3, tri_pts);
        auto is_invalid = checkMinJacDet<15>(nodes_det);
        printf("cand %d eev_col %d, vc = %d is invalid %d\n", cand, eev_col, vc, is_invalid);

        max_invalid = max2(max_invalid, is_invalid);
      }
      invalidities[cand * 2 + eev_col] = max_invalid;
    }
  };
  parallel_for(ncands, f, "coarsen_invalidities");
  auto out_invalid = LOs(invalidities);
  return mesh->sync_subset_array(EDGE, out_invalid, cands2edges, -1, 2);
}

LOs coarsen_invalidities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return coarsen_invalidities_tmpl<3, 3>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return coarsen_invalidities_tmpl<2, 2>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return coarsen_invalidities_tmpl<3, 1>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return coarsen_invalidities_tmpl<2, 1>(mesh, cands2edges, cand_codes);
  }
  OMEGA_H_NORETURN(LOs());
}

Read<I8> filter_coarsen_invalids(
    Read<I8> cand_codes, LOs cand_invalids, LO is_invalid) {
  auto keep_dirs = each_eq_to(cand_invalids, is_invalid);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

}  // end namespace Omega_h
