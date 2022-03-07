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
  auto const old_coords = mesh->coords();

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


  auto v2v_old = mesh->ask_star(0);
  auto v2vv_old = v2v_old.a2ab;
  auto vv2v_old = v2v_old.ab2b;
  auto v2v_new = new_mesh->ask_star(0);
  auto v2vv_new = v2v_new.a2ab;
  auto vv2v_new = v2v_new.ab2b;
  auto const old_v2e = mesh->ask_up(0, 1);
  auto const old_v2ve = old_v2e.a2ab;
  auto const old_ve2e = old_v2e.ab2b;
  
  /*
  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   new_mesh->nverts());
  auto nv2ov_ab2b = new_verts2old_verts.ab2b;
  auto nv2ov_a2ab = new_verts2old_verts.a2ab;
  */
  
  auto find_bdry_edges = OMEGA_H_LAMBDA(LO i) {
    LO v_key = keys2verts[i];
    LO v_onto = keys2verts_onto[i];
    if ((oldvert_gdim[v_key] <= (dim-1)) && (oldvert_gdim[v_onto] <= (dim-1))) {
    //if ((oldvert_gdim[v_key] <= (dim-1)) && (oldvert_gdim[v_onto] == oldvert_gdim[v_key]) && 
      //  (oldvert_gid[v_key] == oldvert_gid[v_onto])) {
      for (LO prod = keys2prods[i]; prod < keys2prods[i+1]; ++prod) {
        LO new_edge = prods2new[prod];
        if ((newedge_gdim[new_edge] <= (dim-1))) {
        //if ((newedge_gdim[new_edge] <= (dim-1)) && (newedge_gdim[new_edge] == oldvert_gdim[v_key])
          //  && (newedge_gid[new_edge] == oldvert_gid[v_key])) {
          auto new_edge_v0 = new_ev2v[new_edge*2 + 0];
          auto new_edge_v1 = new_ev2v[new_edge*2 + 1];
          LO new_edge_old_edge0 = -1;
          LO new_edge_old_edge1 = -1;
          for (LO ve = old_v2ve[v_key]; ve < old_v2ve[v_key+1]; ++ve) {
            auto old_edge = old_ve2e[ve];
            auto old_edge_v0 = old_ev2v[old_edge*2 + 0];
            auto old_edge_v1 = old_ev2v[old_edge*2 + 1];
            LO use_v_old = -1;
            assert((old_edge_v0 == v_key) || (old_edge_v1 == v_key)); 
            if (old_edge_v1 == v_key) {
              use_v_old = old_edge_v0;
            }
            else {
              assert (old_edge_v0 == v_key);
              use_v_old = old_edge_v1;
            }

            if (old_verts2new_verts[use_v_old] == new_edge_v0) {
              new_edge_old_edge0 = old_edge;
            }
            else if (old_verts2new_verts[use_v_old] == new_edge_v1) {
              new_edge_old_edge1 = old_edge;
            }
            else {}
          }
          auto new_old_e0_v0 = old_ev2v[new_edge_old_edge0*2 + 0];
          auto new_old_e0_v1 = old_ev2v[new_edge_old_edge0*2 + 1];
          auto new_old_e1_v0 = old_ev2v[new_edge_old_edge1*2 + 0];
          auto new_old_e1_v1 = old_ev2v[new_edge_old_edge1*2 + 1];
          printf("prob edge %d\n", new_edge_old_edge1);
          printf("newEv {%f %f},{%f %f} made from old {%f %f}, {%f %f}, {%f %f}, {%f %f}\n",
             new_coords[new_edge_v0*dim + 0],new_coords[new_edge_v0*dim + 1], 
             new_coords[new_edge_v1*dim + 0],new_coords[new_edge_v1*dim + 1], 
             old_coords[new_old_e0_v0*dim + 0],old_coords[new_old_e0_v0*dim + 1], 
             old_coords[new_old_e0_v1*dim + 0],old_coords[new_old_e0_v1*dim + 1], 
             old_coords[new_old_e1_v0*dim + 0],old_coords[new_old_e1_v0*dim + 1], 
             old_coords[new_old_e1_v1*dim + 0],old_coords[new_old_e1_v1*dim + 1]); 
          auto c0 = get_vector<dim>(Reals(vert_ctrlPts), new_edge_v0);
          auto c3 = get_vector<dim>(Reals(vert_ctrlPts), new_edge_v1);
          Vector<dim> c1;
          Vector<dim> c2;
          Vector<dim> p1;
          Vector<dim> p2;
          //TODO use 2 pts on 2 collapsing for better curvature
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
          printf("calc midpt %f\n", xi_1);
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
  parallel_for(keys2verts.size(), std::move(find_bdry_edges));

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim, Reals(edge_ctrlPts));
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));

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
         */
    }
  };
  parallel_for(prods2new.size(), std::move(face_blends), "face_blends");
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
