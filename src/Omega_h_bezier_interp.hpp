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
  Few<Real, 6> c1_c2;

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

  c1_c2[0] = Cx[0];
  c1_c2[1] = Cy[0];
  c1_c2[2] = Cz[0];
  c1_c2[3] = Cx[1];
  c1_c2[4] = Cy[1];
  c1_c2[5] = Cz[1];
  return c1_c2;
}

OMEGA_H_DEVICE Few<Real, 9> curve_split_e0_3d(
    Real const new_xi_start, Reals old_vertCtrlPts, Reals old_edgeCtrlPts,
    Int const order, LO const v0_old, LO const v1_old, LO const old_edge) {

  LO const dim = 3;
  OMEGA_H_CHECK(order == 3);
  auto const n_edge_pts = n_internal_ctrlPts(1, order);
  auto c0 = get_vector<3>(old_vertCtrlPts, v0_old);
  auto c1 = vector_3(old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0],
                     old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1],
                     old_edgeCtrlPts[old_edge*n_edge_pts*dim + 2]);
  auto c2 = vector_3(
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0],
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1],
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 2]);
  auto c3 = get_vector<3>(old_vertCtrlPts, v1_old);

  Real const new_xi_3 = new_xi_start + 0.5;
  auto new_c3 = curve_paramToParent_3d(new_xi_3, c0, c1, c2, c3, order);

  Real const new_xi_2 = new_xi_start + xi_2_cube()/2.0;
  auto new_p2 = curve_paramToParent_3d(new_xi_2, c0, c1, c2, c3, order);
  Real const new_xi_1 = new_xi_start + xi_1_cube()/2.0;
  auto new_p1 = curve_paramToParent_3d(new_xi_1, c0, c1, c2, c3, order);

  auto c1_c2 = curve_interpToCtrl_pts_3d(order, c0, new_c3, new_p1, new_p2);
  Few<Real, 9> c1_c2_c3;
  for (LO k = 0; k < dim*2; ++k) {
    c1_c2_c3[k] = c1_c2[k];
  }
  for (LO k = 0; k < dim; ++k) {
    c1_c2_c3[dim*2 + k] = new_c3[k];
  }
  return c1_c2_c3;
}

OMEGA_H_DEVICE Few<Real, 6> curve_split_e1_3d(
    Real const new_xi_start, Reals old_vertCtrlPts, Reals old_edgeCtrlPts,
    Int const order, LO const v0_old, LO const v1_old, LO const old_edge,
    Vector<3> new_c0) {

  LO const dim = 3;
  OMEGA_H_CHECK(order == 3);
  auto const n_edge_pts = n_internal_ctrlPts(1, order);
  auto c0 = get_vector<3>(old_vertCtrlPts, v0_old);
  auto c1 = vector_3(old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0],
                     old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1],
                     old_edgeCtrlPts[old_edge*n_edge_pts*dim + 2]);
  auto c2 = vector_3(
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0],
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1],
      old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 2]);
  auto c3 = get_vector<3>(old_vertCtrlPts, v1_old);

  Real const new_xi_2 = new_xi_start + xi_2_cube()/2.0;
  auto new_p2 = curve_paramToParent_3d(new_xi_2, c0, c1, c2, c3, order);
  Real const new_xi_1 = new_xi_start + xi_1_cube()/2.0;
  auto new_p1 = curve_paramToParent_3d(new_xi_1, c0, c1, c2, c3, order);

  auto c1_c2 = curve_interpToCtrl_pts_3d(order, new_c0, c3, new_p1, new_p2);
  return c1_c2;
}

} // namespace Omega_h

#endif
