#ifndef OMEGA_H_BEZIER_INTERP_HPP
#define OMEGA_H_BEZIER_INTERP_HPP
//TODO rename to bezCurve.hpp and put 2d and 3d fns for curve in here
#include "Omega_h_beziers.hpp"

namespace Omega_h {

OMEGA_H_INLINE Vector<3> curve_paramToParent_3d(
    Real const xi, Vector<3> c0, Vector<3> c1, Vector<3> c2, Vector<3> c3,
    Int const order) {
  Vector<3> new_p;
  for (Int k = 0; k < 3; ++k) {
    new_p[k] = c0[k]*Bi(order, 0, xi) + c1[k]*Bi(order, 1, xi) + 
               c2[k]*Bi(order, 2, xi) + c3[k]*Bi(order, 3, xi);
  }
  return new_p;
}

OMEGA_H_INLINE Few<Real, 6> curve_interpToCtrl_pts_3d(
    Int const order, Vector<3> c0, Vector<3> c3, Vector<3> p1, Vector<3> p2) {
  Few<Real, 6> p1_p2;

  Real const xi_1 = xi_1_cube();
  Real const xi_2 = xi_2_cube();
  auto fx = vector_2(p1[0], p2[0]);
  auto fy = vector_2(p1[1], p2[1]);
  auto fz = vector_2(p1[2], p2[2]);
  auto M1_inv = matrix_2x2(Bi(order, 1, xi_1), Bi(order, 2, xi_1),
                           Bi(order, 1, xi_2), Bi(order, 2, xi_2));
  auto M2 = matrix_2x2(Bi(order, 0, xi_1), Bi(order, 3, xi_1),
                       Bi(order, 0, xi_2), Bi(order, 3, xi_2));
  auto M1 = invert(M1_inv);
  auto cx = vector_2(c0[0], c3[0]);
  auto cy = vector_2(c0[1], c3[1]);
  auto cz = vector_2(c0[2], c3[2]);
  auto Cx = M1*fx - M1*M2*cx;
  auto Cy = M1*fy - M1*M2*cy;
  auto Cz = M1*fz - M1*M2*cz;

  p1_p2[0] = Cx[0];
  p1_p2[1] = Cy[0];
  p1_p2[2] = Cz[0];
  p1_p2[3] = Cx[1];
  p1_p2[4] = Cy[1];
  p1_p2[5] = Cz[1];
  return p1_p2;
}
// workaround CUDA compiler bug
#ifdef OMEGA_H_USE_CUDA
__host__
#endif

#define OMEGA_H_EXPL_INST_DECL(T)
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

} // namespace Omega_h

#endif
