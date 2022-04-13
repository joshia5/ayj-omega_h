#ifndef OMEGA_H_CURVE_VALIDITY_3D_HPP
#define OMEGA_H_CURVE_VALIDITY_3D_HPP

#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

OMEGA_H_DEVICE LO computeTetNodeIndex (LO P, LO i, LO j, LO k) {
  LO l = P-i-j-k;
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

OMEGA_H_DEVICE LO getTetNodeIndex (LO P, LO i, LO j, LO k) {
  return computeTetNodeIndex(P, i, j, k);
}

OMEGA_H_DEVICE Real getTetPartialJacobianDet(Few<Real, 60> nodes, LO P, LO i1,
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

OMEGA_H_DEVICE Real Nijkl(Few<Real, 60> nodes, LO d, LO I, LO J, LO K) {

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
OMEGA_H_DEVICE Few<Real, n> getTetJacDetNodes(LO P, Few<Real, 60> const& elemNodes) {
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

OMEGA_H_INLINE LO checkMinJacDet_3d(Few<Real, 84> const& nodes, LO order) {
  // first 4 vertices
  Real minAcceptable = 0.0;
  for (LO i = 0; i < 4; ++i) {
    if (nodes[i] < minAcceptable) {
      return i+2;
    }
  }

  for (int edge = 0; edge < 6; ++edge){
    for (int i = 0; i < 3*(order-1)-1; ++i){
      if (nodes[4+edge*(3*(order-1)-1)+i] < minAcceptable){
        return 8+edge;
      }
    }
  }

  for (int face = 0; face < 4; ++face){
    for (int i = 0; i < (3*order-4)*(3*order-5)/2; ++i){
      if (nodes[18*order-20+face*(3*order-4)*(3*order-5)/2+i] < minAcceptable) {
          return 14+face;
      }
    }
  }

  for (int i = 0; i < (3*order-4)*(3*order-5)*(3*order-6)/6; ++i) {
    if (nodes[18*order*order-36*order+20+i] < minAcceptable) {
      return 20;
    }
  }

  return -1;
}

LOs checkValidity_3d(Mesh *mesh, LOs new_tets);

} //namespace

#endif
