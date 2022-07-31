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

#include <iostream>
#include <fstream>
namespace Omega_h {

static Few<Real, 10> const BlendedTriangleGetValues(
    Vector<3> const xi, LO b) {
  Few<Real, 10> values;
  Real const blendingTol = 1.e-12;
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  for (LO i = 0; i < 3; ++i)
    values[i] = -std::pow(xii[i], b);
  LO const n = 10;//10 nodes per tri
  for (LO i = 3; i < n; ++i)
    values[i] = 0.0;


  const LO nE = 2;//2 nodes per edge
  Few<LO, 6> tev;
  tev[0] = 0;
  tev[1] = 1;
  tev[2] = 1;
  tev[3] = 2;
  tev[4] = 2;
  tev[5] = 0;

  for (LO i = 0; i < 3; ++i) {
    Real x, xiix, xv;
    Few<Real, 4> v;
    x = xii[tev[i*2 + 0]] + xii[tev[i*2 + 1]];

    if (x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i*2 +1]]/x;

    xv = xiix;
    v[0] = Bi(3, 0, xv);
    v[1] = Bi(3, 3, xv);
    v[2] = Bi(3, 1, xv);
    v[3] = Bi(3, 2, xv);

    for (LO j = 0; j < 2; ++j)
      values[tev[i*2 + j]] += v[j]*std::pow(x, b);
    for (LO j = 0; j < nE; ++j)
      values[3+i*nE+j] = v[2+j]*std::pow(x, b);
  }
  return values;
}

template <Int dim>
void coarsen_curved_verts_and_edges(Mesh *mesh, Mesh *new_mesh, const LOs old2new,
    const LOs prods2new, const LOs old_verts2new_verts, const LOs keys2verts, const LOs keys2verts_onto,
    const LOs keys2prods, const LOs same_verts2new_verts, const LOs same_verts2old_verts) {
  printf("67 coming in fn\n");
  auto const nold_verts = mesh->nverts();
  auto const nold_edges = mesh->nedges();
  auto const nold_faces = mesh->nfaces();
  auto const old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto const old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto const old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto const old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto const old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto const old_v2vf = mesh->ask_up(0,2).a2ab;
  auto const old_vf2f = mesh->ask_up(0,2).ab2b;
  printf("79 ok\n");
  auto const nkeys = keys2verts.size();
  if (!mesh->has_tag(0, "bezier_pts")) {
    mesh->add_tag<Real>(0, "bezier_pts", dim, mesh->coords());
  }
  auto const old_vertCtrlPts = mesh->get_ctrlPts(0);
  auto const old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const old_coords = mesh->coords();
  printf("91 ok\n");

  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_edge = new_mesh->nedges();
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto const nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim, INT8_MAX);
  Write<I8> edge_crv2bdry_dim(nnew_edge, -1);
  Write<I8> edge_dualCone(nnew_edge, -1);
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
    auto const v0 = new_ev2v[e*2 + 0];
    auto const v1 = new_ev2v[e*2 + 1];
    for (LO j=0; j<dim; ++j) {
      edge_ctrlPts[e*n_edge_pts*dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*0.35;
      edge_ctrlPts[e*n_edge_pts*dim + dim + j] = new_coords[v0*dim + j] +
          (new_coords[v1*dim + j] - new_coords[v0*dim + j])*0.65;
    }
  };
  parallel_for(prods2new.size(), std::move(prod_edge_points),
      "prod_edge_points");

  auto const newedge_gdim = new_mesh->get_array<I8>(1, "class_dim");
  auto const newedge_gid = new_mesh->get_array<LO>(1, "class_id");
  auto const oldvert_gdim = mesh->get_array<I8>(0, "class_dim");
  auto const oldvert_gid = mesh->get_array<LO>(0, "class_id");
  auto const oldface_gdim = mesh->get_array<I8>(2, "class_dim");
  auto const oldedge_gdim = mesh->get_array<I8>(1, "class_dim");
  auto const oldface_gid = mesh->get_array<LO>(2, "class_id");

  auto const v2v_old = mesh->ask_star(0);
  auto const v2vv_old = v2v_old.a2ab;
  auto const vv2v_old = v2v_old.ab2b;
  auto const v2v_new = new_mesh->ask_star(0);
  auto const v2vv_new = v2v_new.a2ab;
  auto const vv2v_new = v2v_new.ab2b;
  auto const old_v2e = mesh->ask_up(0, 1);
  auto const old_v2ve = old_v2e.a2ab;
  auto const old_ve2e = old_v2e.ab2b;
  auto vert_ctrlPts_r = Reals(vert_ctrlPts);

  auto new_verts2same_verts = invert_map_by_atomics(same_verts2new_verts,
						    nnew_verts);
  auto const a2ab = new_verts2same_verts.a2ab;
  auto const ab2b = new_verts2same_verts.ab2b;
  Write<LO> nedge_shared_gface_w(nnew_edge, 0);
  Write<Real> cav_edge_len_w(nnew_edge, 0.0);

  if (dim == 3) {
    auto calc_gface_prods = OMEGA_H_LAMBDA(LO i) {
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO const new_edge = prods2new[prod];
        if (newedge_gdim[new_edge] == 2) {
          nedge_shared_gface_w[new_edge] += 1;
          for (LO prod2 = keys2prods[i]; prod2 < keys2prods[i+1]; ++prod2) {
            LO const other_edge = prods2new[prod2];
            if ((other_edge != new_edge) && (newedge_gdim[other_edge] == 2) &&
                (newedge_gid[new_edge] == newedge_gid[other_edge])) {
              //calc new edge length
              LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
              LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
              auto const v0 = get_vector<dim>(new_coords, new_edge_v0);
              auto const v1 = get_vector<dim>(new_coords, new_edge_v1);
              Real length = 0.0;
              for (LO d=0; d<dim; ++d) {
                length += (v1[d] - v0[d])*(v1[d] - v0[d]);
              }
              length = std::sqrt(length);
              cav_edge_len_w[new_edge] += length;
              nedge_shared_gface_w[new_edge] += 1;
            }
          }
          cav_edge_len_w[new_edge] = cav_edge_len_w[new_edge]/
            (nedge_shared_gface_w[new_edge]*1.0);
        }
      }
    };
    parallel_for(keys2verts.size(), std::move(calc_gface_prods));
  }
  auto nedge_shared_gface = Read<LO>(nedge_shared_gface_w);
  auto cav_edge_len = Read<Real>(cav_edge_len_w);

  // for every edge, calc and store 2 unit tangent vectors from either vertex
  Write<Real> tangents(nold_edges*2*dim, 0.0);
  auto calc_tangents = OMEGA_H_LAMBDA (LO e) {
    LO const e_v0 = old_ev2v[e*2 + 0];
    LO const e_v1 = old_ev2v[e*2 + 1];
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

  Write<LO> count_bdry_cavities(1, 0);
  auto count_bdry_cavs = OMEGA_H_LAMBDA(LO i) {
    LO const v_key = keys2verts[i];
    LO const v_onto = keys2verts_onto[i];
    LO has_bdry_cav = -1;
    for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
      LO new_edge = prods2new[prod];
      //edges classified on g_edges or g_faces
      if (((oldvert_gdim[v_key] <= 1) && (oldvert_gdim[v_onto] <= 1) &&
          (newedge_gdim[new_edge] == 1)) || 
         ((dim == 3) && (newedge_gdim[new_edge] == 2))) {
          has_bdry_cav = 1;
        }
    }
    if (has_bdry_cav == 1) 
      atomic_increment(&count_bdry_cavities[0]);
  };
  parallel_for(keys2verts.size(), std::move(count_bdry_cavs));

  Write<I8> edge_cands(nnew_edge, -1);

  auto curve_bdry_edges = OMEGA_H_LAMBDA(LO i) {
    printf("curve bdry pfor i %d\n", i);
    LO const v_key = keys2verts[i];
    LO const v_onto = keys2verts_onto[i];
    Few<I8, 256> bdry_cands_avail;
    Few<Real, 256> all_cands_dist;//max 16 cand and 16 prods
    for (LO k=0; k<256; ++k) all_cands_dist[k] = DBL_MAX;
    for (LO a=0; a<256; ++a) bdry_cands_avail[a] = 1;

    for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
      LO const new_edge = prods2new[prod];
      LO const nedge_shared_gface_i = nedge_shared_gface[new_edge];
      //Real const cav_edge_len_i = cav_edge_len[new_edge];
      LO const new_edge_v0 = new_ev2v[new_edge*2 + 0];
      LO const new_edge_v1 = new_ev2v[new_edge*2 + 1];
      LO const new_edge_v0_old = same_verts2old_verts[ab2b[a2ab[new_edge_v0]]];
      LO const new_edge_v1_old = same_verts2old_verts[ab2b[a2ab[new_edge_v1]]];
      auto const c0 = get_vector<dim>(vert_ctrlPts_r, new_edge_v0);
      auto const c0_coord = get_vector<dim>(old_coords, new_edge_v0_old);
      auto const c3 = get_vector<dim>(vert_ctrlPts_r, new_edge_v1);
      auto const c3_coord = get_vector<dim>(old_coords, new_edge_v1_old);
      auto new_length = 
        (c3_coord[0] - c0_coord[0])*(c3_coord[0] - c0_coord[0]) + 
        (c3_coord[1] - c0_coord[1])*(c3_coord[1] - c0_coord[1]) + 
        (c3_coord[2] - c0_coord[2])*(c3_coord[2] - c0_coord[2]);
      new_length = std::sqrt(new_length);
      Vector<dim> c1,c2;

      //check for new edge, first vertex is vlower or vupper
      LO v_onto_is_first = -1;
      if (v_onto == new_edge_v0_old) {
        v_onto_is_first = 1;
      }
      else {
        OMEGA_H_CHECK (v_onto == new_edge_v1_old);
      }

      //edges classified on g_edges
      if ((oldvert_gdim[v_key] <= 1) && (oldvert_gdim[v_onto] <= 1) &&
          (newedge_gdim[new_edge] == 1)) {
	edge_crv2bdry_dim[new_edge] = 1;
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
	xi_1 = sum_dist1/(sum_dist1 + sum_dist2);
        
	Vector<dim> old_c1;
	for (Int j = 0; j < dim; ++j) {
	  old_c1[j] = (old_p1[j] - B0_quad(xi_1)*c0[j] - 
                       B2_quad(xi_1)*c3[j])/B1_quad(xi_1);
	}
	for (LO d = 0; d < dim; ++d) {
	  c1[d] = (1.0/3.0)*c0[d] + (2.0/3.0)*old_c1[d];
	  c2[d] = (2.0/3.0)*old_c1[d] + (1.0/3.0)*c3[d];
          edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c1[d];
          edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c2[d];
	}
      }
      //end:edges classified on g_edges

      //edges classified on g_faces
      if ((dim == 3) && (newedge_gdim[new_edge] == 2)) {
	edge_crv2bdry_dim[new_edge] = 2;
        //init with straight sided
        for (LO d = 0; d < dim; ++d) {
          c1[d] = edge_ctrlPts[new_edge*n_edge_pts*dim + d];
          c2[d] = edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d];
	}

        //find lower vtx, here upper is v_onto and lower is other end of edge
        Vector<dim> c_lower;
        LO v_lower = -1;
        if (new_edge_v0_old == v_onto) v_lower = new_edge_v1_old;
        if (new_edge_v1_old == v_onto) v_lower = new_edge_v0_old;
        LO e_lower = -1;
        for (LO ve = old_v2ve[v_key]; ve < old_v2ve[v_key + 1]; ++ve) {
          LO const e = old_ve2e[ve];
          if ((v_lower == old_ev2v[e*2+0]) || 
              (v_lower == old_ev2v[e*2+1])) e_lower = e;
        }
        OMEGA_H_CHECK(e_lower);
               
        //count upper edges
        Few<LO, 2> upper_edges;
        Few<I8, 2> from_first_vtx;
        LO count_upper_edge = 0;
        for (LO vf = old_v2vf[v_key]; vf < old_v2vf[v_key + 1]; ++vf) {
          //adj tris of vkey
          LO const adj_t = old_vf2f[vf];
          for (LO te = 0; te < 3; ++te) {
            LO adj_t_e = old_fe2e[adj_t*3 + te];
            //adj edges of tri
            LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
            LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
            //adj verts of edge
            if ((adj_t_e_v0 == v_onto) && (adj_t_e_v1 != v_key)) {
              LO is_duplicate = -1;
              for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = 1;
                  ++count_upper_edge;
                }
              }
            }
            if ((adj_t_e_v1 == v_onto) && (adj_t_e_v0 != v_key)) {
              LO is_duplicate = -1;
              for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = -1;
                  ++count_upper_edge;
                }
              }
            }
            if (oldvert_gdim[v_key] == 1) {
              if ((adj_t_e_v0 == v_onto) && (adj_t_e_v1 == v_key)) {
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = 1;
                  ++count_upper_edge;
                }
              }
              if ((adj_t_e_v1 == v_onto) && (adj_t_e_v0 == v_key)) {
                LO is_duplicate = -1;
                for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
                  if (adj_t_e == upper_edges[upper_e]) is_duplicate = 1;
                }
                if (is_duplicate == -1) {
                  OMEGA_H_CHECK(count_upper_edge < 2);
                  upper_edges[count_upper_edge] = adj_t_e;
                  from_first_vtx[count_upper_edge] = -1;
                  ++count_upper_edge;
                }
              }
            }
          }
        }
        printf("key %d nedgeShareGface %d bdry count_upper_e %d first %d second %d\n", i,
            nedge_shared_gface_i,
            count_upper_edge, upper_edges[0], upper_edges[1]);

        Vector<dim> c_upper;
        Vector<dim> c_avg_upper;
        Few<Real, 2*dim> upper_tangents;
        Vector<dim> t_avg_upper;
        for (LO d = 0; d < dim; ++d) t_avg_upper[d] = 0.0; 
        for (LO upper_e = 0; upper_e < count_upper_edge; ++upper_e) {
          if (from_first_vtx[upper_e] == 1) {
            for (LO d = 0; d < dim; ++d) {
              t_avg_upper[d] += tangents[upper_edges[upper_e]*2*dim + d];
              upper_tangents[upper_e*dim + d] = tangents[upper_edges[upper_e]*2*dim + d];
            }
          }
          if (from_first_vtx[upper_e] == -1) {
            for (LO d = 0; d < dim; ++d) {
              t_avg_upper[d] += tangents[upper_edges[upper_e]*2*dim + dim + d];
              upper_tangents[upper_e*dim + d] = tangents[upper_edges[upper_e]*2*dim + dim + d];
            }
          }
        }

        //calc c using upper avg dir to check concavity
        Real length_t = 0.0;
        LO concave_upper = -1;
        for (LO d = 0; d < dim; ++d) 
          t_avg_upper[d] = t_avg_upper[d]/count_upper_edge;
        for (LO d = 0; d < dim; ++d) 
          length_t += t_avg_upper[d]*t_avg_upper[d]; 
        for (LO d = 0; d < dim; ++d) 
          t_avg_upper[d] = t_avg_upper[d]/std::sqrt(length_t);
        for (LO d = 0; d < dim; ++d) {
          c_avg_upper[d] = old_coords[v_onto*dim + d] + 
                           t_avg_upper[d]*new_length/3.0;
        }
        //calc upper concavity
        Real dist_to_lower = 0.0;
        for (LO d = 0; d < dim; ++d) {
          dist_to_lower += std::pow(
              (c_avg_upper[d] - old_coords[v_lower*dim + d]), 2);
        }
        dist_to_lower = std::sqrt(dist_to_lower);
        if (dist_to_lower > new_length) {
          concave_upper = 1;
          printf("dist %f, leng %f, 1 newE concavity found\n",
              dist_to_lower, new_length);
          for (LO d = 0; d < dim; ++d) {
            c_avg_upper[d] = 
              old_coords[v_onto*dim + d] - t_avg_upper[d]*new_length/3.0;
          }
        }
        //

        if (nedge_shared_gface_i == 1) {
          for (LO d = 0; d < dim; ++d) c_upper[d] = c_avg_upper[d];
        }
        else {
          assert(nedge_shared_gface_i > 1);
          assert(nedge_shared_gface_i <= 32);
          Few<Real, 3*32> cand_tangents;
          Few<Real, 3*32> cand_c;

          Few<Real, 32> cand_angle_to_uppere0;
          Few<LO, 32> sorted_cands;
          for (LO cand = 0; cand < 32; ++cand) {
            sorted_cands[cand] = -1;
            cand_angle_to_uppere0[cand] = DBL_MAX;
          }
          LO const uppere0 = upper_edges[0];
          Vector<dim> uppere0_pt;
          LO const uppere0_v0 = old_ev2v[uppere0*2 + 0];
          LO const uppere0_v1 = old_ev2v[uppere0*2 + 1];
          LO uppere0_vlower = -1;
          if (uppere0_v0 == v_onto) uppere0_vlower = uppere0_v1;
          if (uppere0_v1 == v_onto) uppere0_vlower = uppere0_v0;
          Real uppere0_len = 0.0;
          for (LO d=0; d<dim; ++d) {
            uppere0_len += std::pow((old_coords[uppere0_v0*dim + d] -
                                     old_coords[uppere0_v1*dim + d]), 2);
          }
          uppere0_len = std::sqrt(uppere0_len);
          printf("upper pt\n");
          for (LO d=0; d<dim; ++d) {
            uppere0_pt[d] = old_coords[v_onto*dim + d] + 
              (old_coords[uppere0_vlower*dim + d] - old_coords[v_onto*dim + d])/2.0;
            printf("%f\n", uppere0_pt[d]);
          }

          //calc n locations of new tangent pts
          //printf("newedge leng %f\n", new_length);
          Real upper_theta = acos(
              upper_tangents[0]*upper_tangents[dim + 0] +
              upper_tangents[1]*upper_tangents[dim + 1] +
              upper_tangents[2]*upper_tangents[dim + 2]);
          printf("upper angle %f degree\n", upper_theta*180/PI);
          Vector<dim> n;
          //cross
          n[0] = (upper_tangents[1]*upper_tangents[dim + 2]) - 
                 (upper_tangents[2]*upper_tangents[dim + 1]);
          n[1] = -(upper_tangents[0]*upper_tangents[dim + 2]) +
                 (upper_tangents[2]*upper_tangents[dim + 0]);
          n[2] = (upper_tangents[0]*upper_tangents[dim + 1]) - 
                 (upper_tangents[1]*upper_tangents[dim + 0]);

          for (LO cand = 0; cand < nedge_shared_gface_i; ++cand) {

            //
            Real theta_c = (cand+1)*upper_theta/(nedge_shared_gface_i + 1);
            printf("theta cand %f degree\n",theta_c*180/PI);
            Vector<3> b;
            b[0] = 0.0; b[1] = cos(theta_c); b[2] = cos(upper_theta-theta_c);
            auto A = tensor_3(
                n[0], n[1], n[2],
                upper_tangents[0], upper_tangents[1], upper_tangents[2],
                upper_tangents[dim+0], upper_tangents[dim+1], upper_tangents[dim+2]
                );

            auto A_inv = invert(A);
            auto X = A_inv*b;
            //hc-calc +/- y axis as tang if A is ill-cond
            if (((std::abs(upper_theta - PI) < EPSILON) && 
                ((std::abs(upper_tangents[0]) - 1.0)) < EPSILON)) {
              //+y
              X[0] = upper_tangents[0]*cos(theta_c) +
                     upper_tangents[2]*sin(theta_c);
              X[1] = upper_tangents[1];
              X[2] = -upper_tangents[0]*sin(theta_c) + 
                      upper_tangents[2]*cos(theta_c);
              printf("rotation, n {%f,%f,%f}\n", n[0],n[1],n[2]);
              //-y
              if (n[1] < 0.0) {
                printf("rotation matrix -y \n");
                X[2] = upper_tangents[0]*sin(theta_c) -
                       upper_tangents[2]*cos(theta_c);
              }
            }
            printf("new tang {%f,%f,%f}\n",X[0],X[1],X[2]);
            //

            length_t = 0.0;
            for (LO d=0; d<dim; ++d) {
              cand_tangents[cand*dim + d] = X[d];
              length_t += std::pow(cand_tangents[cand*dim + d], 2);
            }
            printf("#468 newE %d cand_tgts {%f,%f,%f} lenT %f \n", new_edge,
                cand_tangents[cand*dim+0],cand_tangents[cand*dim+1],
                cand_tangents[cand*dim+2],
                length_t);
            for (LO d = 0; d < dim; ++d) {
              cand_tangents[cand*dim + d] = cand_tangents[cand*dim + d]/
                std::sqrt(length_t);
            }
            for (LO d = 0; d < dim; ++d) {
              cand_c[cand*dim + d] = old_coords[v_onto*dim + d] +
                cand_tangents[cand*dim + d]*new_length/3.0;//onto is upper
            }
            if (concave_upper == 1) {
              for (LO d = 0; d < dim; ++d) {
                cand_c[cand*dim + d] = old_coords[v_onto*dim + d] -
                  cand_tangents[cand*dim + d]*new_length/3.0;//onto is upper
              }
              cand_angle_to_uppere0[cand] = 0.0;
              for (LO d = 0; d < dim; ++d) {
                //TODO
                cand_angle_to_uppere0[cand] = acos(
                    upper_tangents[0]*upper_tangents[dim + 0] +
                    upper_tangents[1]*upper_tangents[dim + 1] +
                    upper_tangents[2]*upper_tangents[dim + 2]);
              }
            }
            //
            //for (LO cand = 0; cand < nedge_shared_gface_i; ++cand) {
              //for (LO cand2 = cand; cand2 < nedge_shared_gface_i; ++cand2) {
                //if (cand_angle_to_uppere0[cand] < cand_angle_to_uppere0[cand2]) {
                //  sorted_cands[cand] = cand;
               // }
                //else {
                  //sorted_cands[cand] = cand2;
                  //swap2(cand_angle_to_uppere0[cand] , cand_angle_to_uppere0[cand2]);
                  //swap2(cand_c[cand*dim + 0], cand_c[cand2*dim + 0]);
                  //swap2(cand_c[cand*dim + 1], cand_c[cand2*dim + 1]);
                  //swap2(cand_c[cand*dim + 2], cand_c[cand2*dim + 2]);
                  //swap2(prod_ids[count_p2] , prod_ids[count_prod2_2]);
                //}
              //}
              //printf("sorted cand %d is %d\n", cand,sorted_cands[cand]);
            //}
            //

              printf("#481 candc {%f,%f,%f} cand_tgts {%f,%f,%f} disttoUpp %f newl %f\n", 
                  cand_c[cand*dim+0],cand_c[cand*dim+1],cand_c[cand*dim+2],
                  cand_tangents[cand*dim+0],cand_tangents[cand*dim+1],cand_tangents[cand*dim+2],
                  cand_angle_to_uppere0[cand], new_length);

            //printf("cand dist to uppere0 %f\n", cand_angle_to_uppere0[cand]);
            printf("cand %d tang {%f,%f,%f} upper1 {%f,%f,%f} upper2 {%f,%f,%f}\n", cand, 
                cand_tangents[cand*dim+0],cand_tangents[cand*dim+1],
                cand_tangents[cand*dim+2],
                upper_tangents[0], upper_tangents[1], upper_tangents[2], 
                upper_tangents[dim+0],upper_tangents[dim+1],upper_tangents[dim+2]);
                
          }

          //find dist of all relevant prods to uppere0
          Few<Real, 32> prod_dist_to_uppere0;
          Few<LO, 32> sorted_prods;
          Few<LO, 32> prod_ids;
          for (LO count_p2 = 0; count_p2 < 32; ++count_p2) {
            sorted_prods[count_p2] = -1;
            prod_ids[count_p2] = -1;
            prod_dist_to_uppere0[count_p2] = DBL_MAX;
          }
          LO count_prod2 = 0;
          //TODO find and sort based on angle-to-upper_e0 instead of dist.
          for (LO prod2 = keys2prods[i]; prod2 < keys2prods[i+1]; ++prod2) {
            LO const other_edge = prods2new[prod2];
            if ((newedge_gdim[other_edge] == 2) &&
                (newedge_gid[new_edge] == newedge_gid[other_edge])) {
              LO const oth_edge_v0 = new_ev2v[other_edge*2 + 0];
              LO const oth_edge_v1 = new_ev2v[other_edge*2 + 1];
              LO const oth_edge_v0_old = same_verts2old_verts[ab2b[a2ab[oth_edge_v0]]];
              LO const oth_edge_v1_old = same_verts2old_verts[ab2b[a2ab[oth_edge_v1]]];
              auto const c0_coord2= get_vector<dim>(old_coords, oth_edge_v0_old);
              auto const c3_coord2= get_vector<dim>(old_coords, oth_edge_v1_old);
              Vector<dim> other_p;

              //check for new edge, first vertex is vlower or vupper
              if (v_onto == oth_edge_v0_old) {
                for (LO d = 0; d < dim; ++d)
                  other_p[d] = c0_coord2[d] + (1.0/3.0)*(c3_coord2[d] - c0_coord2[d]);
              }
              else {
                OMEGA_H_CHECK(v_onto == oth_edge_v1_old);
                for (LO d = 0; d < dim; ++d)
                  other_p[d] = c3_coord2[d] + (1.0/3.0)*(c0_coord2[d] - c3_coord2[d]);
              }
              prod_dist_to_uppere0[count_prod2] = 
                std::pow(uppere0_pt[0]-other_p[0], 2) +
                std::pow(uppere0_pt[1]-other_p[1], 2) +
                std::pow(uppere0_pt[2]-other_p[2], 2);
              prod_dist_to_uppere0[count_prod2] = std::sqrt(prod_dist_to_uppere0[count_prod2]);
              prod_ids[count_prod2] = other_edge;
              printf("prod %d dist to uppere0 %f\n", other_edge, prod_dist_to_uppere0[count_prod2]);

              ++count_prod2;
            }
          }
          OMEGA_H_CHECK(count_prod2 == nedge_shared_gface_i);
          //sort prods by dist to upper e0
          for (LO count_p2 = 0; count_p2 < nedge_shared_gface_i; ++count_p2) {
            for (LO count_prod2_2 = count_p2; count_prod2_2 < nedge_shared_gface_i; ++count_prod2_2) {
              if (prod_dist_to_uppere0[count_p2] < prod_dist_to_uppere0[count_prod2_2]) {
                sorted_prods[count_p2] = prod_ids[count_p2];
              }
              else {
                sorted_prods[count_p2] = prod_ids[count_prod2_2];
                swap2(prod_dist_to_uppere0[count_p2] , prod_dist_to_uppere0[count_prod2_2]);
                swap2(prod_ids[count_p2] , prod_ids[count_prod2_2]);
              }
            }
            printf("sorted prod %d is %d\n", count_p2, sorted_prods[count_p2]);
          }
          LO opt_cand_id = -1;
          for (LO count_p2 = 0; count_p2 < nedge_shared_gface_i; ++count_p2) {
            if (prod_ids[count_p2] == new_edge) opt_cand_id = count_p2;
          }
          for (LO d = 0; d < dim; ++d) {
            c_upper[d] = cand_c[opt_cand_id*dim + d];
            //c_upper[d] = cand_c[sorted_cands[opt_cand_id]*dim + d];
          }
          printf("#559 newE 1413 cu {%f,%f,%f} opt_cand_id %d candc {%f,%f,%f}\n", 
              c_upper[0], c_upper[1], c_upper[2],
              opt_cand_id,
              cand_c[opt_cand_id*dim+0],cand_c[opt_cand_id*dim+1],cand_c[opt_cand_id*dim+2]
              );

        }

        Few<LO, 2> lower_edges;
        Few<I8, 2> from_first_vtx_l;
        LO count_lower_edge = 0;
        for (LO vf = old_v2vf[v_key]; vf < old_v2vf[v_key + 1]; ++vf) {
          //adj tris of vkey
          LO const adj_t = old_vf2f[vf];
          for (LO te = 0; te < 3; ++te) {
            LO adj_t_e = old_fe2e[adj_t*3 + te];
            //adj edges of tri
            LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
            LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
            //adj verts of edge
            if ((adj_t_e_v0 == v_lower) && (adj_t_e_v1 != v_key)) {
              LO is_duplicate = -1;
              for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_lower_edge < 2);
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_l[count_lower_edge] = 1;
                  //printf("lower edge n %d is %d\n", count_lower_edge, adj_t_e);
                  ++count_lower_edge;
                }
              }
            }
            if ((adj_t_e_v1 == v_lower) && (adj_t_e_v0 != v_key)) {
              LO is_duplicate = -1;
              for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
                if (adj_t_e == lower_edges[lower_e]) is_duplicate = 1;
              }
              if (is_duplicate == -1) {
                if ((oldface_gid[adj_t] == newedge_gid[new_edge]) && (oldface_gdim[adj_t] == 2)) {
                  OMEGA_H_CHECK(count_lower_edge < 2);
                  lower_edges[count_lower_edge] = adj_t_e;
                  from_first_vtx_l[count_lower_edge] = -1;
                  //printf("lower edge n %d is %d\n", count_lower_edge, adj_t_e);
                  ++count_lower_edge;
                }
              }
            }
          }
        }
        //printf("bdry count_lower_e %d first %d second %d\n", 
        //  count_lower_edge, lower_edges[0], lower_edges[1]);

        Vector<dim> t_lower;
        Few<Real, 2*dim> lower_tangents;
        for (LO d = 0; d < dim; ++d) t_lower[d] = 0.0; 
        for (LO lower_e = 0; lower_e < count_lower_edge; ++lower_e) {
          if (from_first_vtx_l[lower_e] == 1) {
            for (LO d = 0; d < dim; ++d) {
              t_lower[d] += tangents[lower_edges[lower_e]*2*dim + d];
              lower_tangents[lower_e*dim + d] = tangents[lower_edges[lower_e]*2*dim + d];
            }
          }
          if (from_first_vtx_l[lower_e] == -1) {
            for (LO d = 0; d < dim; ++d) {
              t_lower[d] += tangents[lower_edges[lower_e]*2*dim + dim + d];
              lower_tangents[lower_e*dim + d] = tangents[lower_edges[lower_e]*2*dim + dim + d];
            }
          }
        }
        Real lower_cosTheta = (
            lower_tangents[0]*lower_tangents[dim + 0] +
            lower_tangents[1]*lower_tangents[dim + 1] +
            lower_tangents[2]*lower_tangents[dim + 2]);
        printf("lower costheta %f\n", lower_cosTheta);

        OMEGA_H_CHECK(count_lower_edge == 2);
        for (LO d = 0; d < dim; ++d) t_lower[d] = t_lower[d]/count_lower_edge;
        length_t = 0.0;

        for (LO d = 0; d < dim; ++d) length_t += t_lower[d]*t_lower[d]; 
        for (LO d = 0; d < dim; ++d) t_lower[d] = t_lower[d]/std::sqrt(length_t);
        for (LO d = 0; d < dim; ++d) {
          c_lower[d] = old_coords[v_lower*dim + d] + t_lower[d]*new_length/3.0;
        }

        //For concave lower angle//
        {
          Real dist_to_upper = 0.0;
          for (LO d = 0; d < dim; ++d) {
            dist_to_upper += (c_lower[d] - old_coords[v_onto*dim + d])*
              (c_lower[d] - old_coords[v_onto*dim + d]);
          }
          dist_to_upper = std::sqrt(dist_to_upper);
          if (dist_to_upper > new_length) {
            printf("dist %f, leng %f, concavity found\n", dist_to_upper, new_length);
            printf("v_onto {%f,%f,%f} c_lower {%f,%f,%f}\n", old_coords[v_onto*dim + 0],
                old_coords[v_onto*dim + 1], old_coords[v_onto*dim + 2],
                c_lower[0], c_lower[1], c_lower[2]);
            for (LO d = 0; d < dim; ++d) {
              c_lower[d] = old_coords[v_lower*dim + d] - t_lower[d]*new_length/3.0;
            }
            printf(" new c_lower {%f,%f,%f}\n", c_lower[0], c_lower[1], c_lower[2]);
          }
        }
        //

        for (LO d = 0; d < dim; ++d) {
          if (v_onto_is_first == 1) {
            //edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = old_edgeCtrlPts[e_lower*2*dim+d];
            //edge_ctrlPts[new_edge*n_edge_pts*dim + d] = 0.5*(edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] + old_vertCtrlPts[v_onto*dim + d]);
              
            edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_lower[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_upper[d];
          }
          else {
            OMEGA_H_CHECK (v_onto_is_first == -1);
            //edge_ctrlPts[new_edge*n_edge_pts*dim + d] = old_edgeCtrlPts[e_lower*2*dim+d];
            //edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] =
              //0.5*(edge_ctrlPts[new_edge*n_edge_pts*dim + d] + old_vertCtrlPts[v_onto*dim + d]);
            edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_lower[d];
            edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_upper[d];
          }
        }
        if (new_edge == 1413) {
          printf("newE 1413 cu {%f,%f,%f} cl {%f,%f,%f}\n", 
              c_upper[0], c_upper[1], c_upper[2],
              c_lower[0], c_lower[1], c_lower[2]);
        }
      }
    }
    printf("finish pfor i %d\n", i);
  };
  parallel_for(keys2verts.size(), curve_bdry_edges);
  printf("finish bdry crv pfor\n");

  Write<LO> count_dualCone_cavities(1, 0);
  Write<LO> count_interior_dualCone_cavities(1, 0);
  Write<LO> count_interior_cavities(1, 0);
  Write<LO> face_crvVis(nold_faces, -1);
  auto count_dualCone_cav = OMEGA_H_LAMBDA(LO i) {
    LO const v_onto = keys2verts_onto[i];
    LO const v_key = keys2verts[i];
    //printf(" vkey %d vonto %d\n",v_key, v_onto);
    for (LO vf = old_v2vf[v_key]; vf < old_v2vf[v_key + 1]; ++vf) {
      LO const f = old_vf2f[vf];
      face_crvVis[f] = 1;
      //if (v_key == 27) face_crvVis[f] = 1;
    }
    if ((oldvert_gdim[v_key] == dim) && (oldvert_gdim[v_onto] == dim)) {
      atomic_increment(&count_interior_cavities[0]);
    }
    if ((keys2prods[i+1] - keys2prods[i]) == 1) {
      atomic_increment(&count_dualCone_cavities[0]);
      if ((oldvert_gdim[v_key] == dim) && (oldvert_gdim[v_onto] == dim)) {
        atomic_increment(&count_interior_dualCone_cavities[0]);
      }
    }
  };
  parallel_for(keys2prods.size()-1, std::move(count_dualCone_cav));
  //printf("total nkeys %d, nkeys in interior %d nkeys with dual cone cavities %d, %d in interior\n", keys2prods.size()-1, 
    //  count_interior_cavities[0], count_dualCone_cavities[0], count_interior_dualCone_cavities[0]);
  mesh->add_tag<LO>(2, "face_crvVis", 1, Read<LO>(face_crvVis));//old mesh

  if (dim == 3) {
    auto v2t = mesh->ask_up(0, 3);
    auto v2vt = v2t.a2ab;
    auto vt2t = v2t.ab2b;
    auto v2e = mesh->ask_up(0, 1);
    auto v2ve = v2e.a2ab;
    auto ve2e = v2e.ab2b;
    auto te2e = mesh->ask_down(3, 1).ab2b;
 
    auto curve_dualCone_cav = OMEGA_H_LAMBDA (LO i) {
      if ((keys2prods[i+1] - keys2prods[i]) == 1) {
        LO const v_onto = keys2verts_onto[i];
        LO const v_key = keys2verts[i];
        if ((oldvert_gdim[v_key] == dim) && (oldvert_gdim[v_onto] == dim)) {
          LO const new_edge = prods2new[keys2prods[i]];
          edge_dualCone[new_edge] = 1;
          auto new_edge_v0 = new_ev2v[new_edge*2 + 0];
          auto new_edge_v1 = new_ev2v[new_edge*2 + 1];
          auto new_edge_v0_c = get_vector<dim>(new_coords, new_edge_v0);
          auto new_edge_v1_c = get_vector<dim>(new_coords, new_edge_v1);

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
              LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
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
          Real length_t = 0.0;
          for (LO d = 0; d < dim; ++d) length_t += t_avg[d]*t_avg[d]; 
          for (LO d = 0; d < dim; ++d) t_avg[d] = t_avg[d]/std::sqrt(length_t);

          //find lower vtx
          LO v_lower = -1;
          for (LO ve = v2ve[v_key]; ve < v2ve[v_key + 1]; ++ve) {
            //adj edges of vkey
            LO adj_e = ve2e[ve];
            LO adj_e_v0 = old_ev2v[adj_e*2 + 0];
            LO adj_e_v1 = old_ev2v[adj_e*2 + 1];
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
              if (count_not_ring == count_upper_edge) {
                v_lower = other_vtx;
                break;
              }
            }
          }

          Few<LO, 128> lower_edges;
          Few<I8, 128> from_first_vtx_low;
          LO count_lower_edge = 0;
          for (LO vt = v2vt[v_key]; vt < v2vt[v_key + 1]; ++vt) {
            //adj tets of vkey
            LO adj_t = vt2t[vt];
            for (LO te = 0; te < 6; ++te) {
              LO adj_t_e = te2e[adj_t*6 + te];
              //adj edges of tet
              LO adj_t_e_v0 = old_ev2v[adj_t_e*2 + 0];
              LO adj_t_e_v1 = old_ev2v[adj_t_e*2 + 1];
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
          Real length_t_l = 0.0;
          for (LO d = 0; d < dim; ++d) length_t_l += t_avg_l[d]*t_avg_l[d]; 
          for (LO d = 0; d < dim; ++d) t_avg_l[d] = t_avg_l[d]/std::sqrt(length_t_l);
          Vector<dim> c_lower;
          Vector<dim> c_upper;
          auto new_length = 
            (new_edge_v1_c[0] - new_edge_v0_c[0])*(new_edge_v1_c[0] - new_edge_v0_c[0]) + 
            (new_edge_v1_c[1] - new_edge_v0_c[1])*(new_edge_v1_c[1] - new_edge_v0_c[1]) + 
            (new_edge_v1_c[2] - new_edge_v0_c[2])*(new_edge_v1_c[2] - new_edge_v0_c[2]);
          new_length = std::sqrt(new_length);
          for (LO d = 0; d < dim; ++d) {
            c_upper[d] = old_coords[v_onto*dim + d] + t_avg[d]*new_length/3.0;
            c_lower[d] = old_coords[v_lower*dim + d] + t_avg_l[d]*new_length/3.0;
          }

          //check for new edge, first vertex is vlower or vupper
          LO v_onto_is_first = -1;
          if ((std::abs(new_edge_v0_c[0] - old_coords[v_onto*dim + 0]) < EPSILON) && 
              (std::abs(new_edge_v0_c[1] - old_coords[v_onto*dim + 1]) < EPSILON) &&
              (std::abs(new_edge_v0_c[2] - old_coords[v_onto*dim + 2]) < EPSILON)) {
            v_onto_is_first = 1;
          }
          else {
            OMEGA_H_CHECK (
              (std::abs(new_edge_v0_c[0] - old_coords[v_lower*dim + 0]) < EPSILON) && 
              (std::abs(new_edge_v0_c[1] - old_coords[v_lower*dim + 1]) < EPSILON) &&
              (std::abs(new_edge_v0_c[2] - old_coords[v_lower*dim + 2]) < EPSILON));
          }
          
          for (LO d = 0; d < dim; ++d) {
            if (v_onto_is_first == 1) {
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_upper[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_lower[d];
            }
            else {
              OMEGA_H_CHECK (v_onto_is_first == -1);
              edge_ctrlPts[new_edge*n_edge_pts*dim + d] = c_lower[d];
              edge_ctrlPts[new_edge*n_edge_pts*dim + dim + d] = c_upper[d];
            }
          }
        }
      }
    };
    parallel_for(nkeys, std::move(curve_dualCone_cav), "curve_dualCone_cav");
  }

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->add_tag<I8>(1, "edge_crv2bdry_dim", 1, Read<I8>(edge_crv2bdry_dim));
  new_mesh->add_tag<I8>(1, "edge_dualCone", 1, Read<I8>(edge_dualCone));

  return;
}

template <Int dim>
void coarsen_curved_faces(Mesh *mesh, Mesh *new_mesh, const LOs old2new,
    const LOs prods2new) {
  auto const new_fv2v = new_mesh->ask_down(2, 0).ab2b;
  auto const new_fe2e = new_mesh->get_adj(2, 1).ab2b;
  auto const new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto const new_coords = new_mesh->coords();
  auto const nnew_faces = new_mesh->nfaces();
  auto const nnew_edges = new_mesh->nedges();
  auto const nold_faces = mesh->nfaces();
  auto const old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
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
    LO const tri = prods2new[i];
    LO const v0 = new_fv2v[tri*3 + 0];
    LO const v1 = new_fv2v[tri*3 + 1];
    LO const v2 = new_fv2v[tri*3 + 2];
    for (LO j = 0; j < dim; ++j) {
      face_ctrlPts[tri*dim + j] = (new_coords[v0*dim + j] +
         new_coords[v1*dim + j] + new_coords[v2*dim + j])/3.0;
    }
  };
  parallel_for(prods2new.size(), std::move(face_centroids), "face_centroids");

  //loop over new edges
  //if new edge dual cone
  //  find adjacent faces in a for loop
  //  pass face ids as input to blended tri fn 
  //  then take first 9 values and mult them with all ctrl pts of face
  //  that gives the interp pt
  //  face interp to ctrl pt
  if (dim == 3) {
    Vector<3> face_xi;
    face_xi[0] = 1.0/3.0;
    face_xi[1] = 1.0/3.0;
    face_xi[2] = 1.0/3.0;
    auto const weights = BlendedTriangleGetValues(face_xi, 2);
    auto edge_dualCone = new_mesh->get_array<I8>(1, "edge_dualCone");
    auto edge_crv2bdry_dim = new_mesh->get_array<I8>(1, "edge_crv2bdry_dim");
    auto const e2et = new_mesh->ask_up(1, 2).a2ab;
    auto const et2t = new_mesh->ask_up(1, 2).ab2b;
    auto const newedge_gdim = new_mesh->get_array<I8>(1, "class_dim");
    auto const newedge_gid = new_mesh->get_array<LO>(1, "class_id");
    auto const newface_gdim = new_mesh->get_array<I8>(2, "class_dim");
    auto const newface_gid = new_mesh->get_array<LO>(2, "class_id");

    Write<LO> face_crvVis(nnew_faces, -1);
    auto face_blends = OMEGA_H_LAMBDA(LO e) {
      if (edge_dualCone[e] == 1) {
        for (LO et = e2et[e]; et < e2et[e + 1]; ++et) {
          LO const tri = et2t[et];
          face_crvVis[tri] = 1;

          auto p11 = face_blend_interp_3d(3, tri, new_ev2v, new_fe2e,
              vertCtrlPts, edgeCtrlPts, new_fv2v, weights);
          auto newface_c11 = face_interpToCtrlPt_3d(3, tri, new_ev2v, new_fe2e,
              vertCtrlPts, edgeCtrlPts, p11, new_fv2v);
          for (LO k = 0; k < 3; ++k) {
            face_ctrlPts[tri*dim + k] = newface_c11[k];
          }

        }
      }
      if (edge_crv2bdry_dim[e] == 2) {
        for (LO et = e2et[e]; et < e2et[e + 1]; ++et) {
          LO const tri = et2t[et];
          if ((newedge_gid[e] == newface_gid[tri]) && (newface_gdim[tri == 2])) {
            auto p11 = face_blend_interp_3d(3, tri, new_ev2v, new_fe2e,
                vertCtrlPts, edgeCtrlPts, new_fv2v, weights);
            auto newface_c11 = face_interpToCtrlPt_3d(3, tri, new_ev2v, new_fe2e,
                vertCtrlPts, edgeCtrlPts, p11, new_fv2v);
            for (LO k = 0; k < 3; ++k) {
              face_ctrlPts[tri*dim + k] = newface_c11[k];
            }
          }
        }
      }
    };
    parallel_for(nnew_edges, std::move(face_blends), "face_blends");
    new_mesh->add_tag<LO>(2, "face_crvVis", 1, Read<LO>(face_crvVis));
  }

  new_mesh->add_tag<Real>(2, "bezier_pts", dim);
  new_mesh->set_tag_for_ctrlPts(2, Reals(face_ctrlPts));

  return;
}

void check_validity_all_tet(Mesh *new_mesh);
void correct_curved_edges(Mesh *new_mesh);
void check_validity_new_curved_edges(Mesh *new_mesh);

#define OMEGA_H_INST(T)                                                       
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h

#endif
