#ifndef OMEGA_H_CURVE_COARSEN_HPP
#define OMEGA_H_CURVE_COARSEN_HPP

#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_atomics.hpp"

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
  auto const nkeys = keys2verts.size();
  if (!mesh->has_tag(0, "bezier_pts")) {
    mesh->add_tag<Real>(0, "bezier_pts", dim, mesh->coords());
  }
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_coords = mesh->coords();

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<I8> edge_crvtoBdryFace(nnew_edge, -1);
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
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_1_cube();
      edge_ctrlPts[e*n_edge_pts*dim + dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*xi_2_cube();
    }
  };
  parallel_for(prods2new.size(), std::move(prod_edge_points),
      "prod_edge_points");

  auto newedge_gdim = new_mesh->get_array<I8>(1, "class_dim");
  auto newedge_gid = new_mesh->get_array<LO>(1, "class_id");
  auto oldvert_gdim = mesh->get_array<I8>(0, "class_dim");
  auto oldvert_gid = mesh->get_array<LO>(0, "class_id");

  auto v2v_old = mesh->ask_star(0);
  auto v2vv_old = v2v_old.a2ab;
  auto vv2v_old = v2v_old.ab2b;
  auto v2v_new = new_mesh->ask_star(0);
  auto v2vv_new = v2v_new.a2ab;
  auto vv2v_new = v2v_new.ab2b;
  auto const old_v2e = mesh->ask_up(0, 1);
  auto const old_v2ve = old_v2e.a2ab;
  auto const old_ve2e = old_v2e.ab2b;
  auto vert_ctrlPts_r = Reals(vert_ctrlPts);
  //auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
    //                                               nnew_verts);

  auto curve_bdry_edges = OMEGA_H_LAMBDA(LO i) {
    LO v_key = keys2verts[i];
    LO v_onto = keys2verts_onto[i];
    if ((oldvert_gdim[v_key] <= (dim-1)) && (oldvert_gdim[v_onto] <= (dim-1))) {
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO new_edge = prods2new[prod];
        if ((newedge_gdim[new_edge] <= (dim-1))) {
          if (newedge_gdim[new_edge] == 2) {
            edge_crvtoBdryFace[new_edge] = 1;
          }
          auto new_edge_v0 = new_ev2v[new_edge*2 + 0];
          auto new_edge_v1 = new_ev2v[new_edge*2 + 1];

          auto c0 = get_vector<dim>(vert_ctrlPts_r, new_edge_v0);
          auto c3 = get_vector<dim>(vert_ctrlPts_r, new_edge_v1);
          Vector<dim> c1;
          Vector<dim> c2;
          //Vector<dim> p1;
          //Vector<dim> p2;
          auto old_p1 = get_vector<dim>(old_vertCtrlPts, v_key);
          Real xi_1 = 0.5;
          Real sum_dist1 = 0.0;
          Real sum_dist2 = 0.0;
          for (Int j = 0; j < dim; ++j) {
            sum_dist1 += (old_p1[j] - c0[j])*(old_p1[j] - c0[j]);
            sum_dist2 += (old_p1[j] - c3[j])*(old_p1[j] - c3[j]);
          }
          sum_dist1 = std::sqrt(sum_dist1/(dim*1.0));
          sum_dist2 = std::sqrt(sum_dist2/(dim*1.0));
          xi_1 = sum_dist1/(sum_dist1+sum_dist2);
          Vector<dim> old_c1;
          for (Int j = 0; j < dim; ++j) {
            old_c1[j] = (old_p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c3[j])/B1_quad(xi_1);
          }
          for (LO d = 0; d < dim; ++d) {
            c1[d] = (1.0/3.0)*c0[d] + (2.0/3.0)*old_c1[d];
            c2[d] = (2.0/3.0)*old_c1[d] + (1.0/3.0)*c3[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c1[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c2[d];
          }
          break;
        }
      }
    }
  };
  //parallel_for(keys2verts.size(), std::move(curve_bdry_edges));

  Write<LO> count_dualCone_cavities(1, 0);
  Write<LO> count_interior_dualCone_cavities(1, 0);
  auto count_dualCone_cav = OMEGA_H_LAMBDA(LO i) {
    //printf("key %d nprods %d\n",
      //  i, keys2prods[i+1] - keys2prods[i]);
    if ((keys2prods[i+1] - keys2prods[i]) == 1) {
      atomic_increment(&count_dualCone_cavities[0]);

      LO const v_onto = keys2verts_onto[i];
      LO const v_key = keys2verts[i];
      if ((oldvert_gdim[v_key] == dim) && (oldvert_gdim[v_onto] == dim)) {
        atomic_increment(&count_interior_dualCone_cavities[0]);
      }
    }
  };
  parallel_for(keys2prods.size()-1, std::move(count_dualCone_cav));
  printf("total nkeys %d, nkeys with dual cone cavities %d, %d in interior\n", keys2prods.size()-1, 
      count_dualCone_cavities[0], count_interior_dualCone_cavities[0]);

  if (dim == 3) {
    auto v2t = mesh->ask_up(0, 3);
    auto v2vt = v2t.a2ab;
    auto vt2t = v2t.ab2b;
    auto v2e = mesh->ask_up(0, 1);
    auto v2ve = v2e.a2ab;
    auto ve2e = v2e.ab2b;
    auto te2e = mesh->ask_down(3, 1).ab2b;
    auto ev2v = mesh->get_adj(1, 0).ab2b;
 
    // for every edge, calc and store 2 tangents from either vertex
    Write<Real> tangents(nold_edges*2*dim, 0);
    auto calc_tangents = OMEGA_H_LAMBDA (LO e) {
      LO e_v0 = ev2v[e*2 + 0];
      LO e_v1 = ev2v[e*2 + 1];
      Real length1 = 0.0;
      Real length2 = 0.0;
      for (LO d = 0; d < dim; ++d) {
        tangents[e*2*dim + d] = old_edgeCtrlPts[e*2*dim + d] - old_vertCtrlPts[e_v0*dim + d];
        tangents[e*2*dim + dim + d] = old_edgeCtrlPts[e*2*dim + dim + d] - old_vertCtrlPts[e_v1*dim + d];
        length1 += tangents[e*2*dim + d] * tangents[e*2*dim + d];
        length2 += tangents[e*2*dim + dim + d] * tangents[e*2*dim + dim + d];
      }
      for (LO d = 0; d < 3; ++d) {
        tangents[e*2*dim + d] = tangents[e*2*dim + d]/std::sqrt(length1);
        tangents[e*2*dim + dim + d] = tangents[e*2*dim + dim + d]/std::sqrt(length2);
      }
    };
    parallel_for(nold_edges, std::move(calc_tangents));

    auto curve_dualCone_cav = OMEGA_H_LAMBDA (LO i) {
      if ((keys2prods[i+1] - keys2prods[i]) == 1) {
        LO const v_onto = keys2verts_onto[i];
        LO const v_key = keys2verts[i];
        if ((oldvert_gdim[v_key] == dim) && (oldvert_gdim[v_onto] == dim)) {
          LO const new_edge = prods2new[keys2prods[i]];
          auto new_edge_v0 = new_ev2v[new_edge*2 + 0];
          auto new_edge_v1 = new_ev2v[new_edge*2 + 1];
          auto new_edge_v0_c = get_vector<dim>(new_coords, new_edge_v0);
          auto new_edge_v1_c = get_vector<dim>(new_coords, new_edge_v1);
          printf("vkey {%f,%f,%f}, v_onto {%f,%f,%f} new_edge_v0 {%f,%f,%f}, v1 {%f,%f,%f}\n",
              old_coords[v_key*dim + 0], old_coords[v_key*dim + 1], old_coords[v_key*dim + 2],
              old_coords[v_onto*dim + 0], old_coords[v_onto*dim + 1], old_coords[v_onto*dim + 2],
              new_edge_v0_c[0], new_edge_v0_c[1],new_edge_v0_c[2],
              new_edge_v1_c[0], new_edge_v1_c[1],new_edge_v1_c[2]);

          //count upper edges
          Few<LO, 128> upper_edges;
          Few<LO, 128> vtx_ring;
          Few<I8, 128> from_first_vtx;
          LO count_upper_edge = 0;
          for (LO vt = v2vt[v_key]; vt < v2vt[v_key + 1]; ++vt) {
            //adj tets of vkey
            LO adj_t = vt2t[vt];
            for (LO te = 0; te < 6; ++te) {
              LO adj_t_e = te2e[adj_t*6 + te];
              //adj edges of tet
              LO adj_t_e_v0 = ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = ev2v[adj_t_e*2 + 1];
              //adj verts of edge
              if ((adj_t_e_v0 == v_onto) && (adj_t_e_v1 != v_key)) {
                OMEGA_H_CHECK(count_upper_edge < 128);
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = 1;
                  vtx_ring[count_upper_edge] = adj_t_e_v1;
                  ++count_upper_edge;
                }
              }
              if ((adj_t_e_v1 == v_onto) && (adj_t_e_v0 != v_key)) {
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = -1;
                  vtx_ring[count_upper_edge] = adj_t_e_v0;
                  ++count_upper_edge;
                }
              }
            }
          }

          Few<Real, dim> t_avg;
          for (LO d = 0; d < dim; ++d) t_avg[d] = 0.0; 
          for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
            if (from_first_vtx[upper_e] == 1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg[d] += tangents[upper_edges[upper_e]*2*dim + d];
              }
            }
            if (from_first_vtx[upper_e] == -1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg[d] += tangents[upper_edges[upper_e]*2*dim + dim + d];
              }
            }
          }
          for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/count_upper_edge;
          printf("vkey {%f,%f,%f}, v_onto {%f,%f,%f} count upper edges %d tavg {%f,%f,%f} Mag %f \n", 
              old_vertCtrlPts[v_key*dim + 0], old_vertCtrlPts[v_key*dim + 1], old_vertCtrlPts[v_key*dim + 2], 
              old_vertCtrlPts[v_onto*dim + 0], old_vertCtrlPts[v_onto*dim + 1], old_vertCtrlPts[v_onto*dim + 2], 
              count_upper_edge,
              t_avg[0], t_avg[1], t_avg[2],
              std::sqrt(t_avg[0]*t_avg[0] + t_avg[1]*t_avg[1] + t_avg[2]*t_avg[2]));

          //find lower vtx
          LO v_lower = -1;
          for (LO ve = v2ve[v_key]; ve < v2ve[v_key + 1]; ++ve) {
            //adj edges of vkey
            LO adj_e = ve2e[ve];
            LO adj_e_v0 = ev2v[adj_e*2 + 0];
            LO adj_e_v1 = ev2v[adj_e*2 + 1];
            LO other_vtx = -1;
            LO count_not_ring = 0;
            //note id of other end of edge
            if (v_key == adj_e_v0) {
              other_vtx = adj_e_v1;
            }
            if (v_key == adj_e_v1) {
              other_vtx = adj_e_v0;
            }
            if (other_vtx != v_onto) {//eliminate v_onto
              for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                if (other_vtx != vtx_ring[upper_e]) ++count_not_ring;
              }
              if (count_not_ring == count_upper_edge) v_lower = other_vtx;
            }
          }
          printf("v_lower {%f,%f,%f}\n",
            old_coords[v_lower*dim + 0], old_coords[v_lower*dim + 1], old_coords[v_lower*dim + 2]);

          Few<LO, 128> lower_edges;
          Few<I8, 128> from_first_vtx_low;
          LO count_lower_edge = 0;
          for (LO vt = v2vt[v_key]; vt < v2vt[v_key + 1]; ++vt) {
            //adj tets of vkey
            LO adj_t = vt2t[vt];
            for (LO te = 0; te < 6; ++te) {
              LO adj_t_e = te2e[adj_t*6 + te];
              //adj edges of tet
              LO adj_t_e_v0 = ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = ev2v[adj_t_e*2 + 1];
              //adj verts of edge
              if ((adj_t_e_v0 == v_lower) && (adj_t_e_v1 != v_key)) {
                OMEGA_H_CHECK(count_lower_edge < 128);
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = 1;
                  ++count_lower_edge;
                }
              }
              if ((adj_t_e_v1 == v_lower) && (adj_t_e_v0 != v_key)) {
                LO is_duplicate = -1;
                for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                  if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_low[count_lower_edge] = -1;
                  ++count_lower_edge;
                }
              }
            }
          }

          Few<Real, dim> t_avg_l;
          for (LO d = 0; d < dim; ++d) t_avg_l[d] = 0.0; 
          for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
            if (from_first_vtx_low[lower_e] == 1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + d];
              }
            }
            if (from_first_vtx_low[lower_e] == -1) {
              for (LO d = 0; d < dim; ++d) {
                t_avg_l[d] += tangents[lower_edges[lower_e]*2*dim + dim + d];
              }
            }
          }
          for (LO d = 0; d < dim; ++d) t_avg_l[d] = t_avg_l[d]/count_lower_edge;
          printf("count lower edges %d tavg {%f,%f,%f} Mag %f \n", 
              count_lower_edge,
              t_avg_l[0], t_avg_l[1], t_avg_l[2],
              std::sqrt(t_avg_l[0]*t_avg_l[0] + t_avg_l[1]*t_avg_l[1] + t_avg_l[2]*t_avg_l[2]));

        }
      }
    };
    parallel_for(nkeys, std::move(curve_dualCone_cav));
  }

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->add_tag<I8>(1, "edge_crvtoBdryFace", 1, Read<I8>(edge_crvtoBdryFace));

  return;
}

template <Int dim>
void coarsen_curved_faces(Mesh *mesh, Mesh *new_mesh, LOs old2new,
    LOs prods2new) {
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_fe2e = new_mesh->get_adj(2, 1).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_faces = new_mesh->nfaces();
  auto const nold_faces = mesh->nfaces();
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  Write<Real> face_ctrlPts(nnew_faces*1*dim, INT8_MAX);

  auto const vertCtrlPts = mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = mesh->get_ctrlPts(1);

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

  auto face_blends = OMEGA_H_LAMBDA(LO i) {
    LO f = prods2new[i];
    auto v0 = new_fv2v[f*3 + 0];
    auto v1 = new_fv2v[f*3 + 1];
    auto v2 = new_fv2v[f*3 + 2];
    for (LO j=0; j<dim; ++j) {
      face_ctrlPts[f*dim + j] = (
         new_coords[v0*dim + j] +
         new_coords[v1*dim + j] + new_coords[v2*dim + j])/3.0;
      
    /*
    auto e0 = new_fe2e[f*3 + 0];
    auto e1 = new_fe2e[f*3 + 1];
    auto e2 = new_fe2e[f*3 + 2];
    for (LO j=0; j<dim; ++j) {
      face_ctrlPts[f*dim + j] = (
         new_coords[v0*dim + j] +
         new_coords[v1*dim + j] + new_coords[v2*dim + j] +
         edgeCtrlPts[e0*2*dim + j] + 
         edgeCtrlPts[e0*2*dim + dim +j] +
         edgeCtrlPts[e1*2*dim + j] + 
         edgeCtrlPts[e1*2*dim + dim +j] +
         edgeCtrlPts[e2*2*dim + j] + 
         edgeCtrlPts[e2*2*dim + dim +j])/9.0;
    }
    */

    }
  };
  parallel_for(prods2new.size(), std::move(face_blends), "face_blends");
  new_mesh->add_tag<Real>(2, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(2, Reals(face_ctrlPts));

  return;
}
void correct_curved_edges(Mesh *new_mesh);

#define OMEGA_H_INST(T)                                                       
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h

#endif
