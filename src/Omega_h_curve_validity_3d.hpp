#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

OMEGA_H_INLINE LO computeTetNodeIndex (LO P, LO i, LO j, LO k) {
  LO l = P-i-j-k;
  int l = P-i-j-k;
  if(i == P) return 0;
  if(j == P) return 1;
  if(k == P) return 2;
  if(l == P) return 3;
  if(k == 0 && l == 0) return 3+j; // 0-1
  if(i == 0 && l == 0) return 3+(P-1)+k; // 1-2
  if(j == 0 && l == 0) return 3+2*(P-1)+i; // 2-0
  if(j == 0 && k == 0) return 3+3*(P-1)+l; // 0-3
  if(i == 0 && k == 0) return 3+4*(P-1)+l; // 1-3
  if(i == 0 && j == 0) return 3+5*(P-1)+l;// 2-3
  if(l == 0) return k*(P-1)-k*(k-1)/2+j+5*P-2; // 0-1-2
  if(k == 0) return l*(P-1)-l*(l-1)/2+j+5*P-2+(P-2)*(P-1)/2; // 0-1-3
  if(i == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1);// 1-2-3
  if(j == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1)*3/2; // 0-2-3
  return i-P-((i-P+1)*(i-P+2)*(i-P+3))/6+l*(P-1-i)-l*(l-1)/2+k+2*P*P+2;
}

OMEGA_H_INLINE LO getTetNodeIndex (LO P, LO i, LO j, LO k) {
  return computeTetNodeIndex(P, i, j, k);
}

OMEGA_H_INLINE Real getTetPartialJacobianDet(Few<Real, 60> nodes, LO P, LO i1,
    LO j1, LO k1, LO i2, LO j2, LO k2, LO i3, LO j3, LO k3) {
  LO p00 = getTetNodeIndex(P,i1+1,j1,k1);
  LO p01 = getTetNodeIndex(P,i1,j1+1,k1);
  LO p10 = getTetNodeIndex(P,i2+1,j2,k2);
  LO p11 = getTetNodeIndex(P,i2,j2,k2+1);
  LO p20 = getTetNodeIndex(P,i3+1,j3,k3);
  LO p21 = getTetNodeIndex(P,i3,j3,k3);
  return cross(get_vector<3>(nodes, p01)-get_vector<3>(nodes, p00),
               get_vector<3>(nodes, p11)-get_vector<3>(nodes, p10))
         *(get_vector<3>(nodes, p21) - get_vector<3>(nodes, p20));
}

OMEGA_H_INLINE Real Nijkl(Few<Real, 60> nodes, LO d, LO I, LO J, LO k) {

  Real sum = 0.;
  LO CD = quadnomial(3*(d-1),I,J,K);

  for(LO k1 = 0; k1 <= K; ++k1){
    LO k2start = max2(0,K-k1-(d-1));
    for (LO k2 = k2start; k2 <= K-k1; ++k2){
      for (LO j1 = 0; j1 <= J; ++j1){
        LO j2start = max2(0,J-j1-(d-1));
        for (LO j2 = j2start; j2 <= J-j1; ++j2){
          LO i1end = min2(I,d-1-j1-k1);
          for (LO i1 = 0; i1 <= i1end; ++i1){
            LO i2start = max2(0,I+J+K-i1-j1-k1-j2-k2-(d-1));
            LO i2end = min2(I-i1,d-1-j2-k2);
            for (LO i2 = i2start; i2 <= i2end; ++i2){
              LO i3 = I-i1-i2;
              LO j3 = J-j1-j2;
              LO k3 = K-k1-k2;
              sum += quadnomial(d-1,i1,j1,k1)*quadnomial(d-1,i2,j2,k2)
                  *quadnomial(d-1,i3,j3,k3)
                  *getTetPartialJacobianDet(nodes,d,i1,j1,k1,i2,j2,k2,i3,j3,k3);
            }
          }
        }
      }
    }
  }
  return sum*d*d*d/CD;
}

template <Int n>
OMEGA_H_INLINE void getTetJacDetNodes(LO P, Few<Real, 60> const& elemNodes) {
  Few<Real, n> nodes;//n=84
  for (int I = 0; I <= 3*(P-1); ++I) {
    for (int J = 0; J <= 3*(P-1)-I; ++J) {
      for (int K = 0; K <= 3*(P-1)-I-J; ++K) {
        OMEGA_H_CHECK(getTetNodeIndex(3*(P-1),I,J,K) < nodes.size());
        nodes[getTetNodeIndex(3*(P-1),I,J,K)] = Nijkl(elemNodes,P,I,J,K);
      }
    }
  }
  return nodes;
}

  /*
template<Int n>
OMEGA_H_INLINE LO checkMinJacDet(Few<Real, n> const& nodes) {
  // first 3 vertices
  Real minAcceptable = 0.0;
  for (LO i = 0; i < 3; ++i) {
    printf("i %d Nijk[i] %f\n", i, nodes[i]);
    if (nodes[i] < minAcceptable) {
      return i+2;
    }
  }

  Real minJ = 0;
  for (LO edge = 0; edge < 3; ++edge) {
    for (LO i = 0; i < 2*(order-1)-1; ++i) {
      if (nodes[3+edge*(2*(order-1)-1)+i] < minAcceptable) {
        minJ = -1e10;
        // there is no point in doing much with edges if we dont have
        // elevation or subdivision or subdivision matrices
        if (minJ < minAcceptable){
          return 8+edge;
        }
      }
    }
  }

  for (LO i = 0; i < (2*order-3)*(2*order-4)/2; ++i) {
    if (nodes[6*(order-1)+i] < minAcceptable) {
      minJ = -1e10;
      if (minJ < minAcceptable) {
        return 14;
      }
    }
  }
  return -1;
}

LOs checkValidity_3d(Mesh *mesh, LOs new_tris, Int const mesh_dim) {
  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order == 3);

  Write<LO> is_invalid(new_tris.size(), -1);
  //LO const ntri_pts = 10;

  auto check_validity = OMEGA_H_LAMBDA (LO n) {
    //TODO change tri_pts size to 30 for 3d
    Few<Real, 20> tri_pts;//ntri_pts*dim=20
    auto tri = new_tris[n];

    //query the tri's down verts's ctrl pts and store
    for (LO j = 0; j < 3; ++j) {
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
    auto nodes_det = getTriJacDetNodes<15, 2>(order, tri_pts);

    is_invalid[n] = checkMinJacDet<15>(nodes_det);
    if (is_invalid[n] > 0) printf("invalid tri %d\n", tri);
  };
  parallel_for(new_tris.size(), std::move(check_validity));

  return LOs(is_invalid);
}
  */

} //namespace
