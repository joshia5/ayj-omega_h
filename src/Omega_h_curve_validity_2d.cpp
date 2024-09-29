#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

OMEGA_H_INLINE Real calcMinJacDet(Few<Real, 200> nodes, const I8 n_nodes) {
  Real minJ = 1e10;
  for (LO i = 0; i < n_nodes; ++i) {
    printf("i %d detJ %1.15f\n", i, nodes[i]);
    minJ = min2(minJ,nodes[i]);
  }
  return minJ;
}

OMEGA_H_INLINE Real calcMaxJacDet(Few<Real, 200> nodes, const I8 n_nodes) {
  Real maxJ = -1e10;
  for (LO i = 0; i < n_nodes; ++i)
    maxJ = max2(maxJ,nodes[i]);
  return maxJ;
}

LOs checkValidity_2d(Mesh *mesh, LOs new_tris, Int const mesh_dim) {
  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(EDGE);
  auto const n_face_pts = mesh->n_internal_ctrlPts(FACE);
  auto order = mesh->get_max_order();
  LO const nnew_tris = new_tris.size();
  auto const f2v_degree = element_degree(OMEGA_H_SIMPLEX, FACE, VERT);
  auto const f2e_degree = element_degree(OMEGA_H_SIMPLEX, FACE, EDGE);
 
  auto Qs = mesh->ask_qualities();

  Write<LO> is_invalid(nnew_tris, -1);
  //LO const ntri_pts = 10;

  auto check_validity = OMEGA_H_LAMBDA (LO n) {
    Few<Real, 400> tri_pts;//ntri_pts(=200 max)*dim(=2)=400
    //Few<Real, 20> tri_pts;//ntri_pts*dim=20

    auto tri = new_tris[n];

    //query the tri's down verts's ctrl pts and store
    for (LO j = 0; j < 3; ++j) {//3 is tri2vert degree
      if (mesh_dim == 2) {
        auto p = get_vector<2>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
      else {
        OMEGA_H_CHECK (mesh_dim == 3);
        auto p = get_vector<3>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
    }

    //query the tri's down edge's ctrl pts and store

    auto v0 = fv2v[tri*3 + 0];
    auto v1 = fv2v[tri*3 + 1];
    auto v2 = fv2v[tri*3 + 2];
    auto e0 = fe2e[tri*3 + 0];
    auto e1 = fe2e[tri*3 + 1];
    auto e2 = fe2e[tri*3 + 2];
    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    auto flip = vector_3(-1, -1, -1);
    if ((e0v0 == v1) && (e0v1 == v0)) {
      flip[0] = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      flip[1] = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      flip[2] = 1;
    }

    /*
    if (order == 3) {
      for (LO j = 0; j < 3; ++j) {
        LO index = f2v_degree; // after storing 3 pts for vtx
        if (flip[j] == -1) {
          for (I8 d = 0; d < mesh_dim; ++d) {
            tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
              edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
            tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
              edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
          }
        }
        else {
          //for flipped edges
          OMEGA_H_CHECK (flip[j] == 1);
          for (I8 d = 0; d < mesh_dim; ++d) {
            tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
              edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
            tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
              edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
          }
        }
      }

      //query the face's ctrl pt and store
      for (I8 d = 0; d < mesh_dim; ++d) {
        LO index = f2v_degree + f2e_degree*n_edge_pts; // after storing pts for edge
        tri_pts[index*mesh_dim + d] = faceCtrlPts[tri*mesh_dim + d];
      }

      //TODO change to template for mesh_dim
      auto nodes_det = getTriJacDetNodes<15, 2>(order, tri_pts);

      is_invalid[n] = checkMinJacDet<15>(nodes_det, order);
    }
    */

    if (order >= 3) {
    //if (order > 3) {
      OMEGA_H_CHECK(order <= 10);
      for (LO j = 0; j < 3; ++j) {//3 edges per tri
        LO index = f2v_degree; // after storing 3 pts for vertices
        if (flip[j] == -1) {
          for (I8 edge_pt = 0; edge_pt < n_edge_pts; ++edge_pt) {
            for (I8 d = 0; d < mesh_dim; ++d) {
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim +
                      edge_pt*mesh_dim + d] =
              edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim +
                          edge_pt*mesh_dim + d];
            }
          }
        }
        else {
          //for flipped edges
          OMEGA_H_CHECK (flip[j] == 1);
          
          for (I8 edge_pt = 0; edge_pt < n_edge_pts; ++edge_pt) {
            for (I8 d = 0; d < mesh_dim; ++d) {
              tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim +
                edge_pt*mesh_dim + d] =
                edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim +
                (n_edge_pts-1 - edge_pt)*mesh_dim + d];
            }
          }
        }
      }

      //query the face's ctrl pt and store
      for (I8 d = 0; d < mesh_dim; ++d) {
        LO index = f2v_degree + f2e_degree*n_edge_pts; // after storing pts for edge

        for (I8 face_pt = 0; face_pt < n_face_pts; ++face_pt) {
          for (I8 d = 0; d < mesh_dim; ++d) {
            tri_pts[index*mesh_dim + 
                    face_pt*mesh_dim + d] =
              faceCtrlPts[tri*n_face_pts*mesh_dim +
                          face_pt*mesh_dim + d];
          }
        }

      }

      //TODO change to template for mesh_dim
      I8 const detJ_order = 2*(order - 1);
      I8 n_nodes_J = n_total_ctrlPts(FACE, detJ_order);
      //auto nodes_det = getTriJacDetNodes<face_tot_pts, 2>(order, tri_pts);
      auto nodes_det = getTriJacDetNodes<200, 2>(order, tri_pts);
      //auto nodes_det = getTriJacDetNodes<28, 2>(order, tri_pts);
      //28 here is total pts of dim*(order-1) bezier polynomial p=4

      auto const minJ = calcMinJacDet(nodes_det, 3); // only verts
      //auto const minJ = calcMinJacDet(nodes_det, n_nodes_J);
      auto const maxJ = calcMaxJacDet(nodes_det, 3);
      //auto const maxJ = calcMaxJacDet(nodes_det, n_nodes_J);
    printf("tri %d minJ %f, maxJ %f\n",n,minJ, maxJ);

      is_invalid[n] = checkMinJacDet<200>(nodes_det, order);
    }
 
    //if (is_invalid[n] > 0) printf("invalid tri %d\n", tri);

  };
  parallel_for(nnew_tris, std::move(check_validity));

  return LOs(is_invalid);
}

Reals askQuality_2d(Mesh *mesh, LOs new_tris, Int const mesh_dim) {
  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order == 3);
  LO const nnew_tris = new_tris.size();
  
  auto Qs = mesh->ask_qualities();

  Write<Real> Q(nnew_tris, -1e-10);
  //LO const ntri_pts = 10;

  auto check_validity = OMEGA_H_LAMBDA (LO n) {
    Few<Real, 400> tri_pts;//ntri_pts*dim=20
    //Few<Real, 20> tri_pts;//ntri_pts*dim=20
    auto tri = new_tris[n];

    //query the tri's down verts's ctrl pts and store
    for (LO j = 0; j < 3; ++j) {//3 is tri2vert degree
      if (mesh_dim == 2) {
        auto p = get_vector<2>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
      else {
        OMEGA_H_CHECK (mesh_dim == 3);
        auto p = get_vector<3>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
    }

    //query the tri's down edge's ctrl pts and store

    auto v0 = fv2v[tri*3 + 0];
    auto v1 = fv2v[tri*3 + 1];
    auto v2 = fv2v[tri*3 + 2];
    auto e0 = fe2e[tri*3 + 0];
    auto e1 = fe2e[tri*3 + 1];
    auto e2 = fe2e[tri*3 + 2];
    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    auto flip = vector_3(-1, -1, -1);
    if ((e0v0 == v1) && (e0v1 == v0)) {
      flip[0] = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      flip[1] = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      flip[2] = 1;
    }

    for (LO j = 0; j < 3; ++j) {
      LO index = 3;
      if (flip[j] == -1) {
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
        }
      }
      else {
        //for flipped edges
        OMEGA_H_CHECK (flip[j] == 1);
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
        }
      }
    }

    //query the face's ctrl pt and store
    for (I8 d = 0; d < mesh_dim; ++d) {
      LO index = 9;
      tri_pts[index*mesh_dim + d] = faceCtrlPts[tri*mesh_dim + d];
    }

    //TODO change to template for mesh_dim
    auto nodes_det = getTriJacDetNodes<200, 2>(order, tri_pts);

    auto const minJ = calcMinJacDet(nodes_det, 3);
    auto const maxJ = calcMaxJacDet(nodes_det, 3);
    printf("tri %d qs %f, minJ %f, maxJ %f, Q %f\n",n,Qs[n], minJ, maxJ, Q[n]);
    /*
    Q[n] = std::pow((minJ/maxJ), 1./2.)*Qs[n];
    printf("Q=%.10f\n",Q[n]);
    if (Q[n] < 0.001) printf("low quality %f for element %d\n", Q[n], n);
    if (minJ < 0.) {
      Q[n] = 0.;
      //printf("tri %d qs %f, minJ %f, maxJ %f, Q %f\n",n,Qs[n], minJ, maxJ, Q[n]);
    }
    */

  };
  parallel_for(nnew_tris, std::move(check_validity));

  return Reals(Q);
}

} //namespace
